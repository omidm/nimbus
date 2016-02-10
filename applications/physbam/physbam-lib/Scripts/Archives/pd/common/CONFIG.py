# pd Configuration
import os.path

# session directory (location independent) never moved
mon_sessions_file="/data/pd/data/mon_sessions.py"
session_directory="/data/pd/data/sessions"
hosts_filename="/data/pd/data/hosts.py"
# SSL 
keys_directory="/usr/local/adm/pd/keys"
ca_certificate_file=os.path.join(keys_directory,"CA.cert")
server_private_key_file=os.path.join(keys_directory,"server.pkey")
server_certificate_file=os.path.join(keys_directory,"server.cert")
client_private_key_file=os.path.join(keys_directory,"client.pkey")
client_certificate_file=os.path.join(keys_directory,"client.cert")
# pdr server
pdmon_server_host="curvature.stanford.edu"
pdmon_server_port=8886
pdmon_session_file="/data/pd/data/pdmon_sessions.pickle"
# pdhosts server
pdhosts_server_host="curvature.stanford.edu"
pdhosts_server_port=8888
# pdsim server
pdsim_server_host="curvature.stanford.edu"
pdsim_server_port=8889
# pdr server
pdr_server_host="curvature.stanford.edu"
pdr_server_port=9000
# new send
pdsend_server_host="curvature.stanford.edu"
pdsend_server_port=8889
# disk server
pddisk_server_port=8887
