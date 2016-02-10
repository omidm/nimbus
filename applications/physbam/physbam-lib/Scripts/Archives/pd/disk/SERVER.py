#!/usr/bin/python
from pd.common import CONFIG
from pd.common import SOCKET
import os
import mutex
import time
import threading
import socket
import re
valid_path_re=re.compile("^[0-9a-zA-Z_\-\./]+$")
valid_path2_re=re.compile("/solver/vol[0-9]/v/[0-9a-zA-Z_\-\./]+$")

class DISK_SERVER:
    def __init__(self):
        self.hosts_filename=CONFIG.hosts_filename
        self.commands=["List","Create","Destroy","Set_Quota"]
        self.mutex=threading.Lock()
        #self.next_claim_id=1
        #self.Read_Host_List()

    # Public
    def List(self):
        items=[]
        fp=os.popen("zfs list -H")
        #fp=open("/tmp/disk.txt")
        while 1:
            line=fp.readline()
            if line=="": break
            # format is volume, in disk, free, refer, mount point
            items.append(line[:-1].split("\t"))
        return items

    def List_Hash(self):
        volumes=self.List()
        volumehash={}
        for index in xrange(len(volumes)):
            volume=volumes[index][0]
            volumehash[volume]=index
        return volumehash

    def Breakup_VFS(self,vfs):
        vfsparts=vfs.split("/")
        volume,childparts=vfsparts[0],vfsparts[1:]   # i.e. volume='vol0' and vfsparts=['users','aselle']
        parentmount_suffix="--".join(childparts[:-1]) # parentmount_suffix='users'
        childmount_suffix="--".join(childparts) # childmount_suffix='users--aselle'
        parentmount_local="/".join(["",volume,parentmount_suffix])
        vfsmount_local="/".join(["",volume,childmount_suffix])
        vfsmount_client="/".join(["","solver",volume,childmount_suffix])
        parentmount_client="/".join(["","solver",volume,parentmount_suffix])

        return vfsmount_client,vfsmount_local,parentmount_client,parentmount_local

    def Verify_Path(self,directory_name):
        if valid_path_re.match(directory_name)==None:
            return "Invalid characters"
        if valid_path2_re.match(directory_name)==None:
            return "Didn't start with /solver/vol?/v"
        if directory_name.find("..") != -1:
            return "Invalid characters (relative path)"
        if directory_name.find("--") != -1:
            return "Cannot have a path with double dashes"
        if not directory_name.startswith("/solver/"):
            return "Request does not start with /solver/"
        return None

    def Create(self,directory_name):
        verify=self.Verify_Path(directory_name)
        if verify: return verify
        
        vfs=directory_name[len("/solver/"):]
        volhash=self.List_Hash()
        if volhash.has_key(vfs):
            return "volume already exists"
        if not volhash.has_key(os.path.dirname(vfs)):
            return "parent path doesn't exist"
        vfsmount_client,vfsmount_local,parentmount_client,parentmount_local=self.Breakup_VFS(vfs)
        zfs_create="zfs create -o mountpoint=%s %s"%(vfsmount_local,vfs)
        if os.system(zfs_create)!=0:
            return "failed to create vfs"
        else:
            symlink_create="ln -s %s %s"%(vfsmount_client,os.path.join(parentmount_local,os.path.basename(vfs)))
            if os.system(symlink_create)!=0:
                return "failed to create symlink"
            perms_cmd="chmod o+rwt %s"%(vfsmount_local)
            if os.system(perms_cmd):
                return "Failed to set perms"
                
        return vfs
        #print "parentmount '%s' vfsmount '%s'"%(parentmount,vfsmount)

    def Destroy(self,directory_name):
        verify=self.Verify_Path(directory_name)
        if verify: return verify

        vfs=directory_name[len("/solver/"):]
        volhash=self.List_Hash()
        if not volhash.has_key(vfs):
            return "trying to delete volume that does not exist"
        for i in volhash.keys():
            if i.startswith(vfs+"/") and i != vfs:
                return "trying to delete vfs that has children"
        vfsmount_client,vfsmount_local,parentmount_client,parentmount_local=self.Breakup_VFS(vfs)
        symlink_destroy="rm %s"%(os.path.join(parentmount_local,os.path.basename(vfs)))
        zfs_destroy="zfs destroy %s"%vfs
        mountpoint_destroy="rmdir %s"%vfsmount_local
        if os.system(symlink_destroy)!=0:
            return "failed to remove symlink"
        elif os.system(zfs_destroy)!=0:
            return "failed to destroy zfs"
        elif os.system(mountpoint_destroy)!=0:
            return "failed to remove mount point"
        return "done"

    def Set_Quota(self,directory_name,quota):
        verify=self.Verify_Path(directory_name)
        if verify: return verify

        quota_re=re.compile("^[0-9]+[MG]$")
        if not quota_re.match(quota):
            return "invalid quota string"

        vfs=directory_name[len("/solver/"):]
        cmd="zfs set quota=%s %s"%(quota,vfs)
        if os.system(cmd)!=0:
            return "failed to set quota"
        return "quota set"
        
if __name__ == "__main__":
    server=DISK_SERVER()
    hostname=os.popen("hostname").readline()[:-1]
    SOCKET.SERVER(socket.gethostbyname(hostname),CONFIG.pddisk_server_port,server,
                  (CONFIG.server_private_key_file,CONFIG.server_certificate_file,CONFIG.ca_certificate_file))
    pass
