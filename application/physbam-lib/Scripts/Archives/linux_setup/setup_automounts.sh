#!/bin/sh

echo "aperture:/usr/graphics  /usr/graphics           nfs    soft,intr     0 2" >> /etc/fstab
mkdir /usr/graphics
mount -a
mv /etc/auto.u /etc/auto.u-old;ln -s /usr/graphics/adm/etc/auto.u /etc/auto.u

echo "/n      /etc/auto.n     nosuid,intr,soft,rsize=8192,wsize=8192" > /etc/auto.master
echo "/net    /etc/auto.net   nosuid,intr,soft,rsize=8192,wsize=8192" >> /etc/auto.master
echo "/u      /etc/auto.u     nosuid,intr,soft,rsize=8192,wsize=8192" >> /etc/auto.master

echo "cd              -fstype=iso9660,ro,nosuid,nodev :/dev/cdrom" > /etc/auto.misc

echo "$HOST                                           :/" > /etc/auto.n
echo "*                                               &:/" >> /etc/auto.n

/etc/init.d/autofs restart
chkconfig autofs on
