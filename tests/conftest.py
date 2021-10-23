#!/usr/bin/env python3

import pytest
import freezerbox

@pytest.fixture(autouse=True)
def ignore_external_freezerbox_db(monkeypatch):
    mock_db = freezerbox.Database({})
    mock_db.name = "WARNING: ACCESSING EXTERNAL DATABASE"
    mock_db_factory = lambda: mock_db

    monkeypatch.setattr(freezerbox, 'load_db', mock_db_factory)
    monkeypatch.setattr(freezerbox.model, 'load_db', mock_db_factory)

