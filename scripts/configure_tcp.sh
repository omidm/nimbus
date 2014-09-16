#!/bin/bash

sudo echo '# These changes are added to improme the the tcp performance. - omidm' >> /etc/sysctl.conf
sudo echo 'net.core.wmem_max = 536870912' >> /etc/sysctl.conf
sudo echo 'net.core.rmem_max = 536870912' >> /etc/sysctl.conf
sudo echo 'net.ipv4.tcp_rmem = 10240 102400 536870912' >> /etc/sysctl.conf
sudo echo 'net.ipv4.tcp_wmem = 10240 102400 536870912' >> /etc/sysctl.conf
sudo echo 'net.ipv4.tcp_slow_start_after_idle = 0' >> /etc/sysctl.conf
sudo sysctl -p

