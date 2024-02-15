import gzip
import tarfile
import urllib
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from concurrent.futures import ProcessPoolExecutor

import pkg_resources
from tqdm import tqdm as ProgressBar

from .utils import (StopInputGroup, TerminateExecution, env, fileMD5,
                    get_traceback)

g_action_map = {}


def _load_actions():
    global g_action_map  # pylint: disable=global-variable-not-assigned
    for _entrypoint in pkg_resources.iter_entry_points(group="sos_actions"):
        # import actions from entry_points
        # Grab the function that is the actual plugin.
        _name = _entrypoint.name
        try:
            _plugin = _entrypoint.load()
            g_action_map[_name] = _plugin
        except Exception as e:
            from .utils import get_logger

            # look for sos version requirement
            get_logger().warning(f"Failed to load script running action {_entrypoint.name}: {e}")


def sos_run_script(action, script, *args, **kwargs):
    if not g_action_map:
        _load_actions()
    try:
        g_action_map[action](script, *args, **kwargs)
    except KeyError as e:
        raise RuntimeError(f'Undefined script running action {action}') from e


def stop_if(expr, msg="", no_output=False):
    """Abort the execution of the current step or loop and yield
    an warning message `msg` if `expr` is False"""
    if expr:
        raise StopInputGroup(msg=msg, keep_output=not no_output)
    return 0


def done_if(expr, msg=""):
    """Assuming that output has already been generated and stop
    executing the rest of the substep"""
    if expr:
        raise StopInputGroup(msg=msg, keep_output=True)
    return 0


def skip_if(expr, msg=""):
    """Skip the current substep and set _output to empty. Output
    will be removed if already generated."""
    if expr:
        raise StopInputGroup(msg=msg, keep_output=False)
    return 0


def fail_if(expr, msg=""):
    """Raise an exception with `msg` if condition `expr` is False"""
    if expr:
        raise TerminateExecution(msg if msg else "error triggered by action fail_if")
    return 0


def warn_if(expr, msg=""):
    """Yield an warning message `msg` if `expr` is False """
    if expr:
        env.logger.warning(msg)
    return 0


def downloadURL(URL, dest, decompress=False, index=None):
    dest = os.path.abspath(os.path.expanduser(dest))
    dest_dir, filename = os.path.split(dest)
    #
    if not os.path.isdir(dest_dir):
        os.makedirs(dest_dir, exist_ok=True)
    if not os.path.isdir(dest_dir):
        raise RuntimeError(f"Failed to create destination directory to download {URL}")
    #
    message = filename
    if len(message) > 30:
        message = message[:10] + "..." + message[-16:]
    #
    dest_tmp = dest + f".tmp_{os.getpid()}"
    term_width = shutil.get_terminal_size((80, 20)).columns
    try:
        env.logger.debug(f"Download {URL} to {dest}")
        sig = file_target(dest)
        if os.path.isfile(dest):
            prog = ProgressBar(
                desc=message,
                disable=env.verbosity <= 1,
                position=index,
                leave=True,
                bar_format="{desc}",
                total=10000000,
            )
            target = file_target(dest)
            if env.config["sig_mode"] == "build":
                prog.set_description(message + ": \033[32m writing signature\033[0m")
                prog.update()
                target.write_sig()
                prog.close()
                return True
            if env.config["sig_mode"] == "ignore":
                prog.set_description(message + ": \033[32m use existing\033[0m")
                prog.update()
                prog.close()
                return True
            if env.config["sig_mode"] in ("default", "skip", "distributed"):
                prog.update()
                if sig.validate():
                    prog.set_description(message + ": \033[32m Validated\033[0m")
                    prog.update()
                    prog.close()
                    return True
                prog.set_description(message + ":\033[91m Signature mismatch\033[0m")
                target.write_sig()
                prog.update()
        #
        prog = ProgressBar(
            desc=message,
            disable=env.verbosity <= 1,
            position=index,
            leave=True,
            bar_format="{desc}",
            total=10000000,
        )
        #
        # Stop using pycurl because of libcurl version compatibility problems
        # that happen so often and difficult to fix. Error message looks like
        #
        # Reason: Incompatible library version: pycurl.cpython-35m-darwin.so
        # requires version 9.0.0 or later, but libcurl.4.dylib provides version 7.0.0
        #
        # with open(dest_tmp, 'wb') as f:
        #    c = pycurl.Curl()
        #    c.setopt(pycurl.URL, str(URL))
        #    c.setopt(pycurl.WRITEFUNCTION, f.write)
        #    c.setopt(pycurl.SSL_VERIFYPEER, False)
        #    c.setopt(pycurl.NOPROGRESS, False)
        #    c.setopt(pycurl.PROGRESSFUNCTION, prog.curlUpdate)
        #    c.perform()
        # if c.getinfo(pycurl.HTTP_CODE) == 404:
        #    prog.set_description(message + ':\033[91m 404 Error {}\033[0m'.format(' '*(term_width - len(message) - 12)))
        #    try:
        #        os.remove(dest_tmp)
        #    except OSError:
        #        pass
        #    return False
        with open(dest_tmp, "wb") as f:
            try:
                u = urllib.request.urlopen(str(URL))
                try:
                    file_size = int(u.getheader("Content-Length"))
                    prog = ProgressBar(total=file_size, desc=message, position=index, leave=False)
                except Exception:
                    file_size = None
                file_size_dl = 0
                block_sz = 8192
                while True:
                    buffer = u.read(block_sz)
                    if not buffer:
                        break
                    file_size_dl += len(buffer)
                    f.write(buffer)
                    prog.update(len(buffer))
            except urllib.error.HTTPError as e:
                prog.set_description(message + f":\033[91m {e.code} Error\033[0m")
                prog.update()
                prog.close()
                try:
                    os.remove(dest_tmp)
                except OSError:
                    pass
                return False
            except Exception as e:
                prog.set_description(message + f":\033[91m {e}\033[0m")
                prog.update()
                prog.close()
                try:
                    os.remove(dest_tmp)
                except OSError:
                    pass
                return False
        #
        if os.path.isfile(dest):
            os.remove(dest)
        os.rename(dest_tmp, dest)
        decompressed = 0
        if decompress:
            if zipfile.is_zipfile(dest):
                prog.set_description(message + ":\033[91m Decompressing\033[0m")
                prog.update()
                prog.close()
                zfile = zipfile.ZipFile(dest)
                zfile.extractall(dest_dir)
                names = zfile.namelist()
                for name in names:
                    if os.path.isdir(os.path.join(dest_dir, name)):
                        continue
                    if not os.path.isfile(os.path.join(dest_dir, name)):
                        return False
                    decompressed += 1
            elif tarfile.is_tarfile(dest):
                prog.set_description(message + ":\033[91m Decompressing\033[0m")
                prog.update()
                prog.close()
                with tarfile.open(dest, "r:*") as tar:
                    tar.extractall(dest_dir)
                    # only extract files
                    files = [x.name for x in tar.getmembers() if x.isfile()]
                    for name in files:
                        if not os.path.isfile(os.path.join(dest_dir, name)):
                            return False
                        decompressed += 1
            elif dest.endswith(".gz"):
                prog.set_description(message + ":\033[91m Decompressing\033[0m")
                prog.update()
                prog.close()
                decomp = dest[:-3]
                with gzip.open(dest, "rb") as fin, open(decomp, "wb") as fout:
                    buffer = fin.read(100000)
                    while buffer:
                        fout.write(buffer)
                        buffer = fin.read(100000)
                decompressed += 1
        decompress_msg = ("" if not decompressed else
                          f' ({decompressed} file{"" if decompressed <= 1 else "s"} decompressed)')
        prog.set_description(
            message +
            f':\033[32m downloaded{decompress_msg} {" "*(term_width - len(message) - 13 - len(decompress_msg))}\033[0m')
        prog.update()
        prog.close()
        # if a md5 file exists
        # if downloaded files contains .md5 signature, use them to validate
        # downloaded files.
        if os.path.isfile(dest + ".md5"):
            prog.set_description(message + ":\033[91m Verifying md5 signature\033[0m")
            prog.update()
            prog.close()
            with open(dest + ".md5") as md5:
                rec_md5 = md5.readline().split()[0].strip()
                obs_md5 = fileMD5(dest, sig_type='full')
                if rec_md5 != obs_md5:
                    prog.set_description(message + ":\033[91m MD5 signature mismatch\033[0m")
                    prog.update()
                    prog.close()
                    env.logger.warning(
                        f"md5 signature mismatch for downloaded file {filename[:-4]} (recorded {rec_md5}, observed {obs_md5})"
                    )
            prog.set_description(message + ":\033[91m MD5 signature verified\033[0m")
            prog.update()
            prog.close()
    except Exception as e:
        if env.verbosity > 2:
            sys.stderr.write(get_traceback())
        env.logger.error(f"Failed to download: {e}")
        return False
    finally:
        # if there is something wrong still remove temporary file
        if os.path.isfile(dest_tmp):
            os.remove(dest_tmp)
    return os.path.isfile(dest)


def download(URLs, dest_dir=".", dest_file=None, decompress=False, max_jobs=5):
    """Download files from specified URL, which should be space, tab or
    newline separated URLs. The files will be downloaded to specified destination.
    Option "dest_dir" specify the destination directory,
    and "dest_file" specify the output filename, which will otherwise be the same
    specified in the URL. If `filename.md5` files are downloaded, they are used to
    validate downloaded `filename`. If "decompress=True", compressed
    files are decompressed. If `max_jobs` is given, a maximum of `max_jobs`
    concurrent download jobs will be used for each domain. This restriction
    applies to domain names and will be applied to multiple download
    instances.
    """
    if env.config["run_mode"] == "dryrun":
        print(f"HINT: download\n{URLs}\n")
        return None
    if isinstance(URLs, str):
        urls = [x.strip() for x in URLs.split() if x.strip()]
    else:
        urls = list(URLs)

    if not urls:
        env.logger.debug(f"No download URL specified: {URLs}")
        return
    #
    if dest_file is not None and len(urls) != 1:
        raise RuntimeError("Only one URL is allowed if a destination file is specified.")
    #
    if dest_file is None:
        filenames = []
        for idx, url in enumerate(urls):
            token = urllib.parse.urlparse(url)
            # if no scheme or netloc, the URL is not acceptable
            if not all([getattr(token, qualifying_attr) for qualifying_attr in ("scheme", "netloc")]):
                raise ValueError(f"Invalid URL {url}")
            filename = os.path.split(token.path)[-1]
            if not filename:
                raise ValueError(f"Cannot determine destination file for {url}")
            filenames.append(os.path.join(dest_dir, filename))
    else:
        token = urllib.parse.urlparse(urls[0])
        if not all([getattr(token, qualifying_attr) for qualifying_attr in ("scheme", "netloc")]):
            raise ValueError(f"Invalid URL {url}")
        filenames = [dest_file]
    #
    succ = [(False, None) for x in urls]
    with ProcessPoolExecutor(max_workers=max_jobs) as executor:
        for idx, (url, filename) in enumerate(zip(urls, filenames)):
            # if there is alot, start download
            succ[idx] = executor.submit(downloadURL, url, filename, decompress, idx)
    succ = [x.result() for x in succ]

    # for su, url in zip(succ, urls):
    #    if not su:
    #        env.logger.warning('Failed to download {}'.format(url))
    failed = [y for x, y in zip(succ, urls) if not x]
    if failed:
        if len(urls) == 1:
            raise RuntimeError("Failed to download {urls[0]}")
        raise RuntimeError(f"Failed to download {failed[0]} ({len(failed)} out of {len(urls)})")
    return 0
