#!/usr/bin/python

import os

session_directory=os.environ["PSIM_SESSION_DIRECTORY"]
session_id=os.environ["PSIM_SESSION_ID"]

os.system("screen -d -R psim_"+str(session_id))
