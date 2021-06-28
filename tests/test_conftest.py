#!/usr/bin/env python3

def test_ignore_external_freezerbox_db():
    import freezerbox
    db = freezerbox.load_db()
    assert db.name == "WARNING: ACCESSING EXTERNAL DATABASE"
