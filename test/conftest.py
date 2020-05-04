import pytest
import uuid
import os
import tempfile


@pytest.fixture
def config_factory():
    filename = tempfile.NamedTemporaryFile(suffix='.yml', delete=False).name

    def get_config(text):
        with open(filename, 'w') as conf:
            conf.write(text)
        return filename

    yield get_config
    #os.remove(filename)
