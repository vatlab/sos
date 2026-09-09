#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import operator

from .targets import BaseTarget
from .utils import env, textMD5

# tested in this order so that two-character operators are matched before the
# single-character operators that they start with
_VERSION_OPERATORS = (
    ("==", operator.eq),
    ("!=", operator.ne),
    ("<=", operator.le),
    (">=", operator.ge),
    ("<", operator.lt),
    (">", operator.gt),
)


class Py_Module(BaseTarget):
    """A target for a Python module."""

    LIB_STATUS_CACHE: dict[str, str] = {}

    def __init__(self, module, version=None, autoinstall=False):
        super().__init__()
        if not isinstance(module, str):
            raise ValueError("A string is expected for module name.")
        self._module = module.strip()
        self._version = version.strip() if isinstance(version, str) else version
        for opt in ("==", ">=", ">", "<=", "<", "!="):
            if opt in self._module:
                if self._version is not None:
                    raise ValueError(f"Specifying 'version=' option in addition to '{module}' is not allowed")
                self._module, self._version = (x.strip() for x in self._module.split(opt, 1))
                if "," in self._version:
                    raise ValueError(f"SoS does not yet support multiple version comparisons. {module} provided")
                self._version = opt + self._version
                break
        self._autoinstall = autoinstall

    def _version_satisfied(self, ver):
        """Check if version ver of an installed module satisfies the requested version,
        which is either a plain version or a version prefixed by a comparison operator
        such as '>=1.0'."""
        from packaging.version import InvalidVersion
        from packaging.version import parse as parse_version

        compare = operator.eq
        requested = self._version
        for opt, opt_compare in _VERSION_OPERATORS:
            if self._version.startswith(opt):
                compare = opt_compare
                requested = self._version[len(opt) :]
                break
        else:
            if self._version[0] in ("=", ">", "<", "!"):
                # an incomplete comparison operator such as '=1.0'
                return False
        try:
            return compare(parse_version(ver), parse_version(requested))
        except InvalidVersion as e:
            env.logger.debug(f"Failed to compare version {ver} against {self._version}: {e}")
            return False

    def _install(self, name, autoinstall):
        """Check existence of Python module and install it using command
        pip install if necessary."""
        import importlib
        from importlib import metadata

        spam_spec = importlib.util.find_spec(name)
        reinstall = False
        if spam_spec is not None and self._version:
            mod = importlib.__import__(name)
            ver = getattr(mod, "__version__", None)
            if ver is None:
                try:
                    ver = metadata.version(name)
                except Exception as e:
                    env.logger.debug(f"Failed to get version of {name}: {e}")
            if ver is None:
                env.logger.warning(f"Cannot determine version of installed {name} to compare against {self._version}.")
                reinstall = True
            else:
                env.logger.debug(f"Comparing existing version {ver} against requested version {self._version}")
                if not self._version_satisfied(ver):
                    env.logger.warning(
                        f"Version {ver} of installed {name} does not match specified version {self._version}."
                    )
                    reinstall = True
        if spam_spec and not reinstall:
            return True
        if not autoinstall:
            return False
        # try to install it?
        import subprocess

        cmd = (
            ["pip", "install"]
            + ([] if self._version else ["-U"])
            + [
                self._module + (self._version if self._version else "")
                if self._autoinstall is True
                else self._autoinstall
            ]
        )
        env.logger.info(f"Installing python module {name} with command {' '.join(cmd)}")
        ret = subprocess.call(cmd)
        if reinstall:
            import sys

            importlib.reload(sys.modules[name])
        # try to check version
        return ret == 0 and self._install(name, False)

    def target_exists(self, mode="any"):
        if (self._module, self._version) in self.LIB_STATUS_CACHE:
            return self.LIB_STATUS_CACHE[(self._module, self._version)]
        ret = self._install(self._module, self._autoinstall)
        self.LIB_STATUS_CACHE = {x: y for x, y in self.LIB_STATUS_CACHE.items() if x[0] != self._module}
        self.LIB_STATUS_CACHE[(self._module, self._version)] = ret
        return ret

    def target_name(self):
        return self._module

    def target_signature(self, mode="any"):
        # we are supposed to get signature of the module, but we cannot
        return textMD5("Python module " + self._module)
