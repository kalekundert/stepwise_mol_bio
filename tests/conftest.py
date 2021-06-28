#!/usr/bin/env python3

import pytest
import freezerbox

@pytest.fixture(autouse=True)
def ignore_external_freezerbox_db(monkeypatch):
    monkeypatch.setattr(
            freezerbox, 'load_db',
            lambda: freezerbox.Database("WARNING: ACCESSING EXTERNAL DATABASE"),
    )

