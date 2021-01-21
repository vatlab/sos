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
    except subprocess.CalledProcessedError as e:
        return (e.returncode, e.output)
    except Exception as e:
        return (1, f"error: failed to check output of {cmd}: {e}")
    finally:
        if orig_dir is not None:
            os.chdir(orig_dir)


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
                action = decode_msg(server_socket.recv())
                env.logger.info(f'RECV: {action}')
                if action == 'alive':
                    reply_msg = "yes"
                elif action[0] == 'signature':
                    reply_msg = handle_signature(*actions[1:])
                elif action[0] == 'exists':
                    reply_msg = handle_exists(*actions[1:])
                elif action[0] == 'check_output':
                    reply_msg = handle_check_output(*actions[1:])
                else:
                    reply_msg = f'Unrecognized request {action}'
                env.logger.info(f'SEND: {str(reply_msg)}[:40]')
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
