# AGENTS.md

This repo uses automation agents (local or CI) to keep code healthy and consistent. Keep this document short and actionable.

## Targets (what agents must ensure)
- **Package management**: use `pipenv` for environment and dependency management, but maintain `pyproject.toml` (preferred) and keep `setup.py` / `setup.cfg` in sync if present.
- **Python**: use **Python 3.13**. If the project is not yet on 3.13, upgrade it and CI accordingly.
- **Static checks**: enforce `mypy` and `flake8` on all tracked Python files.
- **Tests**: run `pytest` with `pytest-cov`; fail if coverage drops below the configured threshold.
- **Docs**: use **Sphinx docstring style**. When you touch a file, add/refresh docstrings.
- **Dependencies**: keep them up to date with minimal, safe upgrades.
- **Changes**: when you touch code, you also add/update tests.
- **Cleanup**: remove unused and unexposed code.


## Environment & packaging
- **Env manager**: use **pipenv** for local workflow and CI execution.
- **Project metadata**: prefer a single **pyproject.toml** as the source of truth.
- **Legacy files**: if `setup.py`/`setup.cfg` exist, keep them minimal and synchronized with `pyproject.toml`, or remove them once migration is complete.
- **Locking**: commit `Pipfile`, `Pipfile.lock` and `requirements.txt`.

## Quickstart (local)
```bash
# Python 3.13 environment
pyenv install -s 3.13
pyenv local 3.13

# Install pipenv
pip install pipenv

# Create and use pipenv environment
pipenv --python 3.13
pipenv install --dev    # dev deps include: mypy, flake8, pytest, pytest-cov

# Static analysis
pipenv run flake8 .
pipenv run mypy .

# Tests + coverage
pipenv run pytest --cov=codfreq --cov-report=term-missing

# Update Pipfile.lock
pipenv lock --dev --clear

# Generate requirements.txt
make requirements.txt
```

## CI expectations (example)
- Use Python 3.13 runner.
- Steps:
  1. `pip install -e .[dev]`
  2. `flake8 .`
  3. `mypy .`
  4. `pytest --cov=codfreq --cov-report=xml` (record artifact, enforce threshold)

## Sphinx docstring style (minimal rules)
- Use Sphinx fields: `:param name:`, `:type name:`, `:returns:`, `:rtype:`, `:raises:`.
- Document public functions, classes, and non-trivial private helpers you modify.
- Keep docstrings accurate to code behavior and types.

## Dependency hygiene
```bash
# With pipenv
pipenv update           # safe minor/patch upgrades per constraints
pipenv update <pkg>

# If using pip/requirements:
pip list --outdated
pip install -U <pkg>
# If using pip-tools:
pip-compile --upgrade
pip-sync
```
- Pin in requirements/lockfile as appropriate.
- Run tests and type checks after any upgrade.

## Pre-commit (recommended)
```bash
pip install pre-commit
pre-commit install
# Example hooks: flake8, mypy (via local hook), trailing-whitespace, end-of-file-fixer
```

## Pull request checklist
- [ ] Code runs on Python 3.13.
- [ ] Added/updated tests for all touched code.
- [ ] `flake8` and `mypy` pass.
- [ ] `pytest --cov` meets coverage threshold.
- [ ] Docstrings updated using Sphinx style.
- [ ] Dependencies updated if relevant and tests still pass.
