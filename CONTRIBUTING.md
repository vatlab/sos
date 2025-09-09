# Contributing to SoS

Thank you for your interest in contributing to SoS! This document provides guidelines and instructions for contributing to the project.

## Table of Contents

- [Development Setup](#development-setup)
- [Building the Project](#building-the-project)
- [Running Tests](#running-tests)
- [Code Style](#code-style)
- [Making Changes](#making-changes)
- [Submitting Pull Requests](#submitting-pull-requests)
- [Deployment](#deployment)

## Development Setup

### Prerequisites

- Python 3.6 or higher (Python 3.8+ recommended)
- Git
- pip or conda

### Setting Up Your Development Environment

1. **Fork and clone the repository**

```bash
# Fork the repository on GitHub, then:
git clone https://github.com/YOUR_USERNAME/SoS.git
cd SoS
git remote add upstream https://github.com/vatlab/SoS.git
```

2. **Create a virtual environment** (recommended)

```bash
# Using venv
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Or using conda
conda create -n sos-dev python=3.9
conda activate sos-dev
```

3. **Install in development mode**

```bash
# Install the package in editable mode with development dependencies
pip install -e ".[dev]"

# For platform-specific dependencies
pip install -e ".[dev,unix]"  # On Linux/macOS
pip install -e ".[dev,win]"   # On Windows

# For all optional dependencies
pip install -e ".[dev,dot,unix]"
```

4. **Install pre-commit hooks** (optional but recommended)

```bash
pip install pre-commit
pre-commit install
```

## Building the Project

SoS now uses modern Python packaging with `pyproject.toml` and the `hatchling` build backend.

### Building Distributions

```bash
# Install build tool
pip install build

# Build both wheel and source distribution
python -m build

# Output will be in dist/
ls dist/
# sos-0.25.2-py3-none-any.whl
# sos-0.25.2.tar.gz
```

### Building Documentation

```bash
# Install documentation dependencies
pip install sphinx sphinx-rtd-theme

# Build docs
cd docs
make html
# View at docs/_build/html/index.html
```

## Running Tests

### Running All Tests

```bash
# Change to test directory
cd test

# Build required Docker containers for Docker tests
sh build_test_docker.sh

# Run all tests
python run_tests.py

# Or using pytest directly
pytest

# Run with coverage
pytest --cov=sos --cov-report=html
```

### Running Specific Tests

```bash
# Run a specific test file
pytest test/test_actions.py

# Run a specific test function
pytest test/test_actions.py::test_function_name

# Run tests matching a pattern
pytest -k "test_bash"

# Run tests with verbose output
pytest -v

# Stop on first failure
pytest -x
```

### Running Linting and Type Checks

```bash
# Run pylint
python -m pylint --rcfile .github/linters/.python-lint src

# Run pre-commit checks
pre-commit run --all-files

# Format code with black (if configured)
black src/

# Sort imports with isort (if configured)
isort src/
```

## Code Style

- Follow PEP 8 guidelines
- Use meaningful variable and function names
- Add docstrings to all public functions and classes
- Keep functions small and focused
- Write tests for new features

### Example Code Style

```python
def calculate_sum(numbers: List[int]) -> int:
    """Calculate the sum of a list of numbers.
    
    Args:
        numbers: List of integers to sum
        
    Returns:
        The sum of all numbers
        
    Raises:
        ValueError: If the list is empty
    """
    if not numbers:
        raise ValueError("Cannot sum an empty list")
    return sum(numbers)
```

## Making Changes

### Workflow

1. **Create a new branch**

```bash
git checkout -b feature/your-feature-name
# or
git checkout -b fix/issue-number
```

2. **Make your changes**

```bash
# Edit files
# Add tests for new functionality
# Update documentation if needed
```

3. **Test your changes**

```bash
# Run relevant tests
pytest test/test_your_changes.py

# Run linting
pre-commit run --all-files
```

4. **Commit your changes**

```bash
git add .
git commit -m "Brief description of changes"
```

### Commit Message Guidelines

- Use present tense ("Add feature" not "Added feature")
- Use imperative mood ("Move cursor to..." not "Moves cursor to...")
- Limit first line to 72 characters
- Reference issues and pull requests when relevant

Example:
```
Fix Docker execution with Python 3.12

- Update subprocess calls for compatibility
- Add error handling for missing containers
- Update tests to cover new behavior

Fixes #1234
```

## Submitting Pull Requests

1. **Push your branch**

```bash
git push origin feature/your-feature-name
```

2. **Create a Pull Request**

- Go to https://github.com/vatlab/SoS
- Click "New Pull Request"
- Select your branch
- Fill in the PR template with:
  - Description of changes
  - Related issues
  - Testing performed
  - Checklist items

3. **PR Requirements**

- [ ] Tests pass locally
- [ ] Code follows project style guidelines
- [ ] Documentation updated if needed
- [ ] Commit messages are clear
- [ ] PR description explains the changes

## Deployment

### Publishing to PyPI

For maintainers only:

1. **Update version**

```bash
# Edit src/sos/_version.py
# Update __version__ = "X.Y.Z"
```

2. **Create a release commit**

```bash
git add src/sos/_version.py
git commit -m "Release version X.Y.Z"
git tag vX.Y.Z
git push origin master --tags
```

3. **Build and upload**

```bash
# Clean previous builds
rm -rf dist/ build/

# Build distributions
python -m build

# Check the distributions
twine check dist/*

# Upload to TestPyPI first (optional)
twine upload --repository testpypi dist/*

# Upload to PyPI
twine upload dist/*
```

### Publishing to Conda-Forge

The conda-forge package is maintained separately. After PyPI release:

1. Fork https://github.com/conda-forge/sos-feedstock
2. Update the version and SHA256 in `recipe/meta.yaml`
3. Submit a PR to conda-forge

## Getting Help

- **Issues**: https://github.com/vatlab/SoS/issues
- **Discussions**: https://github.com/vatlab/SoS/discussions
- **Gitter Chat**: https://gitter.im/vatlab/SoS
- **Documentation**: https://vatlab.github.io/sos-docs

## Additional Resources

- [SoS Documentation](https://vatlab.github.io/sos-docs)
- [Extending SoS](https://vatlab.github.io/sos-docs/doc/user_guide/extending_sos.html)
- [Python Packaging Guide](https://packaging.python.org)
- [Modern Python Packaging](https://packaging.python.org/en/latest/tutorials/packaging-projects/)

## License

By contributing to SoS, you agree that your contributions will be licensed under the same 3-clause BSD License that covers the project.