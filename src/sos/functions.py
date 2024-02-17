from .utils import (StopInputGroup, TerminateExecution, env)

# g_action_map = {}

# def _load_actions():
#     global g_action_map  # pylint: disable=global-variable-not-assigned
#     for _entrypoint in pkg_resources.iter_entry_points(group="sos_actions"):
#         # import actions from entry_points
#         # Grab the function that is the actual plugin.
#         _name = _entrypoint.name
#         try:
#             _plugin = _entrypoint.load()
#             g_action_map[_name] = _plugin
#         except Exception as e:
#             from .utils import get_logger

#             # look for sos version requirement
#             get_logger().warning(f"Failed to load script running action {_entrypoint.name}: {e}")

# def sos_run_script(action, script, *args, **kwargs):
#     '''Call script-execution actions.'''
#     if not g_action_map:
#         _load_actions()
#     try:
#         g_action_map[action](script, *args, **kwargs)
#     except KeyError as e:
#         raise RuntimeError(f'Undefined script running action {action}') from e


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
