#!/usr/bin/python
import sys
import os
import dialog
import CLIENT_LIBRARY


d=dialog.Dialog(dialog="dialog")
d.add_persistent_args(["--backtitle", "PSIM Status"])

while 1:
    id=CLIENT_LIBRARY.Get_Session(d)
    CLIENT_LIBRARY.Status(d,id)
