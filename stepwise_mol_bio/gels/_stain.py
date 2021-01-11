#!/usr/bin/env python3

import autoprop
from inform import Error, did_you_mean
from stepwise_mol_bio import Main

@autoprop
class Stain(Main):

    def get_protocol(self):
        p = self.staining_protocol
        if self.image_type:
            p += self.imaging_protocol
        return p

    def get_staining_protocol(self):
        try:
            return self.staining_protocols[self.stain_type]()
        except KeyError:
            raise Error(f"unknown staining protocol {self.stain_type!r}, did you mean {did_you_mean(self.stain_type, self.staining_protocols)!r}")

    def get_staining_protocols(self):
        raise NotImplementedError

    def get_imaging_protocol(self):
        key = self.image_type
        try:
            return self.imaging_protocols[key]()
        except KeyError:
            raise Error(f"unknown imaging protocol {key!r}, did you mean {did_you_mean(key, self.imaging_protocols.keys())!r}")

    def get_imaging_protocols(self):
        raise NotImplementedError

