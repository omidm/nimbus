from pd.common import SOCKET
from pd.common import CONFIG

def host_client():
    print "HOST CLIENT IS OBSOLETE"
    sys.exit(1)
    return SOCKET.CLIENT(CONFIG.pdhosts_server_host,CONFIG.pdhosts_server_port,
              (CONFIG.client_private_key_file,CONFIG.client_certificate_file,CONFIG.ca_certificate_file))

def send_client(timeout=0):
    return SOCKET.CLIENT(CONFIG.pdsend_server_host,CONFIG.pdsend_server_port,None,timeout)
    #return SOCKET.CLIENT(CONFIG.pdsend_server_host,CONFIG.pdsend_server_port,(CONFIG.client_private_key_file,CONFIG.client_certificate_file,CONFIG.ca_certificate_file),timeout)
             # (CONFIG.client_private_key_file,CONFIG.client_certificate_file,CONFIG.ca_certificate_file))

def sim_client():
    print "SIM CLIENT IS OBSOLETE"
    sys.exit(1)
    return SOCKET.CLIENT(CONFIG.pdsim_server_host,CONFIG.pdsim_server_port,
              (CONFIG.client_private_key_file,CONFIG.client_certificate_file,CONFIG.ca_certificate_file))

def disk_client(volume):
    print "DISK CLIENT IS CURRENTLY OFFLINE"
    sys.exit(1)

    volume_servers={"vol0":"solverh1","vol1":"solverh1","vol2":"solverh2","vol3":"solverh2"}
    return SOCKET.CLIENT(volume_servers[volume],CONFIG.pddisk_server_port,
              (CONFIG.client_private_key_file,CONFIG.client_certificate_file,CONFIG.ca_certificate_file))
    
def mon_client(timeout=0):
    return SOCKET.CLIENT(CONFIG.pdmon_server_host,CONFIG.pdmon_server_port,None,timeout)
    #return SOCKET.CLIENT(CONFIG.pdmon_server_host,CONFIG.pdmon_server_port,(CONFIG.client_private_key_file,CONFIG.client_certificate_file,CONFIG.ca_certificate_file),timeout)
#    return SOCKET.CLIENT(CONFIG.pdmon_server_host,CONFIG.pdmon_server_port,
#              (CONFIG.client_private_key_file,CONFIG.client_certificate_file,CONFIG.ca_certificate_file),timeout)



