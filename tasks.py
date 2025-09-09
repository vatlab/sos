"""
Invoke tasks for SoS development.

Usage:
    invoke --list       # List all available tasks
    invoke test         # Run tests
    invoke format       # Format code
    invoke lint         # Run linting
    invoke clean        # Clean build artifacts
"""

from pathlib import Path

from invoke import task

# Project paths
PROJECT_ROOT = Path(__file__).parent
SRC_DIR = PROJECT_ROOT / "src"
TEST_DIR = PROJECT_ROOT / "test"
DOCS_DIR = PROJECT_ROOT / "docs"


@task
def format(c, check=False):
    """Format code using ruff."""
    cmd = "ruff format"
    if check:
        cmd += " --check"
    cmd += f" {SRC_DIR} {TEST_DIR} tasks.py"
    
    print(f"{'Checking' if check else 'Formatting'} code with ruff...")
    c.run(cmd, pty=True)


@task
def lint(c, fix=False):
    """Run linting with ruff."""
    cmd = "ruff check"
    if fix:
        cmd += " --fix"
    cmd += f" {SRC_DIR} {TEST_DIR} tasks.py"
    
    print("Running ruff linter...")
    c.run(cmd, pty=True)


@task
def test(c, verbose=False, coverage=False, markers="", keyword="", failfast=False):
    """
    Run tests with pytest.
    
    Args:
        verbose: Show verbose test output
        coverage: Run with coverage reporting
        markers: Run tests matching given markers (e.g., "not slow")
        keyword: Run tests matching given keyword expression
        failfast: Stop on first failure
    """
    print("Running tests...")
    
    with c.cd(TEST_DIR):
        # Build test docker images if needed
        if (TEST_DIR / "build_test_docker.sh").exists():
            print("Building test Docker containers...")
            c.run("sh build_test_docker.sh", pty=True)
        
        # Build pytest command
        cmd = "pytest"
        
        if verbose:
            cmd += " -v"
        
        if coverage:
            cmd += " --cov=sos --cov-report=html --cov-report=term"
        
        if markers:
            cmd += f" -m '{markers}'"
        
        if keyword:
            cmd += f" -k '{keyword}'"
        
        if failfast:
            cmd += " -x"
        
        # Run tests
        c.run(cmd, pty=True)
        
        if coverage:
            print("\nCoverage report generated in htmlcov/index.html")


@task
def test_file(c, file, verbose=False):
    """Run tests for a specific file."""
    print(f"Running tests for {file}...")
    cmd = f"pytest {file}"
    if verbose:
        cmd += " -v"
    c.run(cmd, pty=True)


@task
def clean(c, all=False):
    """
    Clean build artifacts and caches.
    
    Args:
        all: Remove all generated files including .venv
    """
    print("Cleaning build artifacts...")
    
    patterns = [
        "build/",
        "dist/",
        "*.egg-info",
        "**/__pycache__",
        "**/*.pyc",
        "**/*.pyo",
        ".pytest_cache/",
        ".coverage",
        "htmlcov/",
        ".ruff_cache/",
        "**/.ipynb_checkpoints",
    ]
    
    if all:
        patterns.extend([".venv/", "uv.lock"])
    
    for pattern in patterns:
        c.run(f"rm -rf {pattern}", warn=True)
    
    print("✓ Cleaned build artifacts")


@task
def build(c):
    """Build distribution packages."""
    print("Building distribution packages...")
    
    # Clean old builds first
    clean(c)
    
    # Build with uv if available, otherwise use build module
    try:
        c.run("which uv", hide=True, warn=True)
        print("Using uv to build...")
        c.run("uv build", pty=True)
    except:
        print("Using python -m build...")
        c.run("python -m build", pty=True)
    
    # List generated files
    c.run("ls -la dist/", pty=True)


@task
def install(c, dev=False, extras=""):
    """
    Install the package.
    
    Args:
        dev: Install in development mode with dev dependencies
        extras: Additional extras to install (comma-separated)
    """
    print("Installing package...")
    
    if dev:
        # Development installation
        try:
            c.run("which uv", hide=True, warn=True)
            print("Using uv for development installation...")
            extras_cmd = "--all-extras" if not extras else f"--extra {extras.replace(',', ' --extra ')}"
            c.run(f"uv sync {extras_cmd}", pty=True)
        except:
            print("Using pip for development installation...")
            extras_str = "[dev]" if not extras else f"[dev,{extras}]"
            c.run(f"pip install -e '.{extras_str}'", pty=True)
    else:
        # Regular installation
        extras_str = f"[{extras}]" if extras else ""
        c.run(f"pip install -e '.{extras_str}'", pty=True)


@task
def docs(c, serve=False, port=8000):
    """
    Build documentation.
    
    Args:
        serve: Start a local server to preview docs
        port: Port for the documentation server
    """
    if not DOCS_DIR.exists():
        print("No docs directory found. Skipping documentation build.")
        return
    
    with c.cd(DOCS_DIR):
        print("Building documentation...")
        c.run("make clean", warn=True)
        c.run("make html", pty=True)
        
        if serve:
            print(f"Serving documentation at http://localhost:{port}")
            with c.cd("_build/html"):
                c.run(f"python -m http.server {port}", pty=True)


@task
def check(c):
    """Run all checks (format check, lint, tests)."""
    print("Running all checks...\n")
    
    print("=" * 60)
    print("Checking code formatting...")
    print("=" * 60)
    format(c, check=True)
    
    print("\n" + "=" * 60)
    print("Running linter...")
    print("=" * 60)
    lint(c)
    
    print("\n" + "=" * 60)
    print("Running tests...")
    print("=" * 60)
    test(c)
    
    print("\n✓ All checks passed!")


@task
def pre_commit(c):
    """Run pre-commit hooks on all files."""
    print("Running pre-commit hooks...")
    c.run("pre-commit run --all-files", pty=True)


@task
def deps_update(c, package=""):
    """
    Update dependencies.
    
    Args:
        package: Specific package to update (updates all if not specified)
    """
    print("Updating dependencies...")
    
    try:
        c.run("which uv", hide=True, warn=True)
        if package:
            print(f"Updating {package}...")
            c.run(f"uv lock --upgrade-package {package}", pty=True)
        else:
            print("Updating all dependencies...")
            c.run("uv lock --upgrade", pty=True)
        c.run("uv sync", pty=True)
    except:
        print("uv not found. Please install uv for dependency management.")
        print("Install with: curl -LsSf https://astral.sh/uv/install.sh | sh")


@task
def deps_show(c, outdated=False):
    """
    Show project dependencies.
    
    Args:
        outdated: Show only outdated packages
    """
    if outdated:
        print("Checking for outdated packages...")
        c.run("uv pip list --outdated", pty=True)
    else:
        print("Current dependencies:")
        c.run("uv pip list", pty=True)


@task
def version(c, bump=""):
    """
    Show or bump version.
    
    Args:
        bump: Version bump type (major, minor, patch)
    """
    version_file = SRC_DIR / "sos" / "_version.py"
    
    if not bump:
        # Just show current version
        with open(version_file) as f:
            for line in f:
                if line.startswith("__version__"):
                    version = line.split("=")[1].strip().strip('"')
                    print(f"Current version: {version}")
                    return
    else:
        # Bump version
        import re
        
        with open(version_file) as f:
            content = f.read()
        
        # Extract current version
        match = re.search(r'__version__ = "(\d+)\.(\d+)\.(\d+)"', content)
        if not match:
            print("Could not parse version from _version.py")
            return
        
        major, minor, patch = map(int, match.groups())
        
        if bump == "major":
            major += 1
            minor = 0
            patch = 0
        elif bump == "minor":
            minor += 1
            patch = 0
        elif bump == "patch":
            patch += 1
        else:
            print(f"Invalid bump type: {bump}. Use major, minor, or patch.")
            return
        
        new_version = f"{major}.{minor}.{patch}"
        
        # Update version in file
        content = re.sub(
            r'__version__ = "[^"]+"',
            f'__version__ = "{new_version}"',
            content
        )
        
        with open(version_file, "w") as f:
            f.write(content)
        
        print(f"Version bumped to {new_version}")
        print("Don't forget to commit and tag the version change!")


@task
def release(c, test_pypi=False):
    """
    Create a release and upload to PyPI.
    
    Args:
        test_pypi: Upload to TestPyPI instead of PyPI
    """
    print("Preparing release...")
    
    # Run all checks first
    check(c)
    
    # Build packages
    build(c)
    
    # Upload to PyPI
    if test_pypi:
        print("Uploading to TestPyPI...")
        c.run("twine upload --repository testpypi dist/*", pty=True)
        print("\nTest installation with:")
        print("  pip install --index-url https://test.pypi.org/simple/ sos")
    else:
        print("Uploading to PyPI...")
        c.run("twine upload dist/*", pty=True)
        print("\n✓ Release uploaded to PyPI!")


@task(name="list")
def list_tasks(c):
    """List all available tasks."""
    c.run("invoke --list", pty=True)


# Aliases for common tasks
@task
def fmt(c, check=False):
    """Alias for format."""
    format(c, check=check)


@task
def l(c, fix=False):
    """Alias for lint."""
    lint(c, fix=fix)


@task
def t(c, verbose=False, coverage=False):
    """Alias for test."""
    test(c, verbose=verbose, coverage=coverage)


@task
def c(c, all=False):
    """Alias for clean."""
    clean(c, all=all)