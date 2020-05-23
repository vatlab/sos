import pytest
import uuid
import os
import tempfile
import pathlib


@pytest.fixture
def config_factory():
    filename = tempfile.NamedTemporaryFile(suffix='.yml', delete=False).name

    def get_config(text):
        with open(filename, 'w') as conf:
            conf.write(text)
        return filename

    yield get_config
    os.remove(filename)


@pytest.fixture
def tempfile_factory():

    temp_files = []

    def get_tempfiles(names):
        if isinstance(names, str):
            names = [names]
        for name in names:
            pathlib.Path(name).touch()
            temp_files.append(name)
        return temp_files

    yield get_tempfiles

    for temp_file in temp_files:
        try:
            os.remove(temp_file)
        except Exception:
            pass
