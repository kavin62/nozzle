from typing import Generator

def walk(self) -> Generator: ...
def body_line_iterator(msg, decode: bool = ...) -> Generator: ...
def typed_subpart_iterator(msg, maintype=..., subtype=...) -> Generator: ...