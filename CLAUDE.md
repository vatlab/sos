# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

SoS (Script of Scripts) is a polyglot workflow system that consists of:
- **SoS Workflow**: A workflow engine for executing workflows in both process- and outcome-oriented styles
- **SoS Notebook**: A Jupyter-based polyglot notebook allowing multiple kernels in one notebook

This repository contains the SoS Workflow engine implementation.

## Development Commands

### Installation
```bash
# Development installation with uv (recommended)
uv venv
source .venv/bin/activate
uv sync --all-extras

# Traditional pip installation
pip install -e ".[dev]"
```

### Using Invoke Tasks (Recommended)
```bash
# Show all available tasks
invoke --list

# Common development tasks
invoke format      # Format code with ruff
invoke lint        # Check code style
invoke test        # Run tests
invoke check       # Run all checks (format, lint, test)
invoke clean       # Clean build artifacts
invoke build       # Build distribution packages

# Test specific scenarios
invoke test --verbose --coverage
invoke test --keyword "test_name"
invoke test --markers "not slow"
invoke test-file test/test_actions.py

# Dependency management
invoke deps-show --outdated
invoke deps-update
invoke deps-update --package pytest

# Version management
invoke version              # Show current version
invoke version --bump patch # Bump version (major/minor/patch)

# Release
invoke release              # Build and upload to PyPI
invoke release --test-pypi  # Upload to TestPyPI
```

### Manual Testing
```bash
# Run all tests
cd test && python run_tests.py

# Run specific test file
pytest test/test_actions.py

# Run specific test function
pytest test/test_actions.py::test_function_name

# Build test Docker containers (required for Docker tests)
cd test && sh build_test_docker.sh
```

### Code Quality
```bash
# Using invoke (recommended)
invoke format && invoke lint && invoke test

# Manual commands
ruff check src/
ruff format src/
pre-commit run --all-files
```

## Code Architecture

### Core Components

**Workflow Engine** (`src/sos/`)
- `__main__.py`: Entry point for `sos` command-line interface
- `workflow_executor.py`: Main workflow execution logic
- `targets.py`: Target types (file_target, sos_step, dynamic, etc.)
- `dag.py`: Directed acyclic graph for workflow dependencies
- `controller.py`: Workflow controller and job management

**Actions** (`src/sos/actions*.py`)
- Language-specific action implementations (bash, python, R, julia, etc.)
- Docker and Singularity container support
- Script execution with various interpreters

**Task System**
- `task_engines.py`: Task execution engines (process, PBS/cluster systems)
- `task_executor.py`: Task execution and management
- Remote task execution support

**Extensions**
- Entry points system for plugins (targets, actions, task engines, previewers)
- Language modules installed separately (sos-r, sos-python, etc.)

### Key Concepts

1. **Steps and Substeps**: Workflows consist of numbered steps (e.g., `[10]`, `[20]`) that can have substeps when processing multiple inputs
2. **Input/Output Groups**: Steps process inputs in groups with parameters like `group_by`, `paired_with`, `for_each`
3. **Targets**: Various target types including files, sos_step dependencies, dynamic targets
4. **Task Distribution**: Support for local and remote execution via task engines

## Important Notes

- Python 3.6+ required
- Uses pytest for testing framework
- Extensive use of entry points for plugin architecture
- Docker and Singularity support for containerized execution
- Cross-platform support (Linux, macOS, Windows)