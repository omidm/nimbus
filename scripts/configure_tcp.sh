#!/bin/bash

sudo echo '# These changes are added to improme the the tcp performance. - omidm'
sudo sysctl -w net.core.wmem_max=536870912
sudo sysctl -w net.core.rmem_max=536870912
sudo sysctl -w net.ipv4.tcp_rmem='10240 102400 536870912'
sudo sysctl -w net.ipv4.tcp_wmem='10240 102400 536870912'
sudo sysctl -w net.ipv4.tcp_slow_start_after_idle=0

