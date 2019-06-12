#setup file mostly copied from official shiny_server git repo
yum clean all -y --enablerepo=*
yum install -y epel-release
yum update -y --disablerepo=epel

# Enable EPEL
rpm -Uvh https://dl.fedoraproject.org/pub/epel/7/x86_64/Packages/e/epel-release-7-11.noarch.rpm 

# On this minimal install, we need wget
yum install -y coreutils
yum install -y yum-utils
yum install -y wget
yum install -y which
yum install -y bzip2
yum install -y git
yum install -y gcc
yum install -y fontconfig
yum install -y libcurl libcurl-devel
yum install -y openssl-devel 

# Install R
yum install R -y

wget https://s3.amazonaws.com/rstudio-shiny-server-os-build/centos-6.3/x86_64/VERSION -O "version.txt"
VERSION=`cat version.txt`

# Install the latest SS build
wget "https://s3.amazonaws.com/rstudio-shiny-server-os-build/centos-6.3/x86_64/shiny-server-$VERSION-rh6-x86_64.rpm" -O ss-latest.rpm
yum install --nogpgcheck ss-latest.rpm -y

sudo su - \
    -c "R -e \"install.packages(c('shiny', 'httpuv', 'rmarkdown', 'devtools', 'RJDBC'), repos='http://cran.rstudio.com/')\""

sudo R -e 'devtools::install_github("tidyverse/ggplot2")'

systemctl disable firewalld 
systemctl stop firewalld
sed -i 's/enforcing/disabled/g' /etc/selinux/config
systemctl enable shiny-server
systemctl start shiny-server
