#!/usr/bin/env python3

import autoprop
from inform import Error, did_you_mean

@autoprop
class Stain:
    default_image_type = None

    def __init__(self):
        self.stain_type = None
        self.image_type = None

    @classmethod
    def main(cls, *args, **kwargs):
        try:
            from docopt import docopt
            args = docopt(*args, **kwargs)
            stain = cls.from_docopt(args)
            print(stain.protocol)
        except Error as err:
            err.report()

    @classmethod
    def from_docopt(cls, args):
        return cls()

    def get_protocol(self):
        return self.staining_protocol + self.imaging_protocol

    def get_staining_protocol(self):
        try:
            return self.staining_protocols[self.stain_type]()
        except KeyError:
            raise Error(f"unknown staining protocol {self.stain_type!r}, did you mean {did_you_mean(self.stain_type, self.staining_protocols)!r}")

    def get_staining_protocols(self):
        raise NotImplementedError

    def get_imaging_protocol(self):
        key = self.image_type or self.default_image_type or ''
        try:
            return self.imaging_protocols[key]()
        except KeyError:
            raise Error(f"unknown imaging protocol {key!r}, did you mean {did_you_mean(key, self.imaging_protocols.keys())!r}")

    def get_imaging_protocols(self):
        raise NotImplementedError

