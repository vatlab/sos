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
# Development installation
pip install -e .
pip install -r requirements_dev.txt

# Full installation with language modules
pip install sos sos-pbs sos-notebook sos-bash sos-python sos-r
```

### Testing
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

### Linting and Code Quality
```bash
# Run pylint
python -m pylint --rcfile .github/linters/.python-lint src

# Run pre-commit hooks
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