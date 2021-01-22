import os
import sys
import subprocess

import zmq
from .messages import encode_msg, decode_msg
from .utils import env
from .targets import sos_targets, file_target

# used by eval
assert sos_targets
assert file_target


def handle_signature(items, target_dir):
    try:
        orig_dir = os.getcwd()
        os.chdir(target_dir)

        items = eval(items)

        return str(sos_targets(items).target_signature())
    except Exception as e:
        return f"error: {e}"
    finally:
        os.chdir(orig_dir)


def handle_exists(items, target_dir):
    try:
        orig_dir = os.getcwd()
        os.chdir(target_dir)

        items = eval(items)
        return "yes" if sos_targets(items).target_exists() else "no"
    except Exception as e:
        return f"error: {e}"
    finally:
        os.chdir(orig_dir)


def handle_check_output(cmd, workdir, kwargs):
    try:
        if workdir:
            orig_dir = os.getcwd()
            os.chdir(workdir)
        else:
            orig_dir = None
        output = subprocess.check_output(cmd, shell=True, **kwargs).decode()
        return (0, output)
    except subprocess.CalledProcessError as e:
        return (e.returncode, e.output)
    except Exception as e:
        return (1, f"error: failed to check output of {cmd}: {e}")
    finally:
        if orig_dir is not None:
            os.chdir(orig_dir)


g_running_procs = {}

def handle_check_call(cmd, workdir, kwargs):
    global g_running_procs
    try:
        if workdir:
            orig_dir = os.getcwd()
            os.chdir(workdir)
        else:
            orig_dir = None

        p = subprocess.Popen(cmd, shell=isinstance(cmd, str), **kwargs)
        g_running_procs[p.pid] = p

        return 'running', p.pid
    except Exception as e:
        return 'exception', e
    finally:
        if orig_dir is not None:
            os.chdir(orig_dir)

def handle_poll_call(pid):
    global g_running_procs
    if pid not in g_running_procs:
        return 'exception', ValueError(f'Invalid call id {pid}')

    ret = g_running_procs[pid].poll()

    if ret is None:
        return 'running', ''

    return 'done', ret

def cmd_server(args, workflow_args):
    if workflow_args:
        raise RuntimeError(f'Unrecognized arguments {" ".join(workflow_args)}')

    context = zmq.Context()
    server_socket = context.socket(zmq.REP)
    server_socket.bind(f"tcp://*:{args.port}")

    env.verbosity = args.verbosity

    try:
        while True:
            #  Wait for next request from client
            if server_socket.poll(
                    -1 if args.duration is None else 1000 * args.duration,
                    zmq.POLLIN):
                params = decode_msg(server_socket.recv())
                env.logger.info(f'RECV: {params}')
                if params == 'alive':
                    reply_msg = "yes"
                elif params[0] == 'signature':
                    reply_msg = handle_signature(*params[1:])
                elif params[0] == 'exists':
                    reply_msg = handle_exists(*params[1:])
                elif params[0] == 'check_output':
                    reply_msg = handle_check_output(*params[1:])
                elif params[0] == 'check_call':
                    reply_msg = handle_check_call(*params[1:])
                elif params[0] == 'poll_call':
                    reply_msg = handle_poll_call(*params[1:])
                else:
                    reply_msg = f'Unrecognized request {params}'
                env.logger.info(f'SEND: {str(reply_msg)[:40]}')
                server_socket.send(encode_msg(reply_msg))
            else:
                break
    except Exception as e:
        env.logger.error(e)
        sys.exit(1)
    finally:
        server_socket.close()
    # after idling args.duration, quit
    sys.exit(0)
