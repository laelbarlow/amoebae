#!/usr/bin/env bash
# Copyright 2018 Lael D. Barlow
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 

# Script for setting up AMOEBAE to run in a singularity container after cloning
# the git repository. This is a bash script so it does not depend on python3
# being installed (for a setup.py script to work) and it does not depend on
# developer tools being installed on MacOS (for make to work).


# Start timing this script.
SECONDS=0


#######################################################
# Copy the example settings file.
if [ ! -f "settings.py" ]; then
    printf "\tWriting settings file.\n\n"
    cp settings.py.example settings.py
fi


#######################################################
# Make scripts executable.
chmod a+x amoebae
chmod a+x misc_scripts/*.py
chmod a+x misc_scripts/*.sh


#######################################################
# If Vagrant is not installed, detect what operating system is being used, and
# install Vagrant if necessary.

if [ "$(uname)" == "Darwin" ]; then
    printf "\nDetected MacOS."
    # Do something under Mac OS X platform

    # Check whether Vagrant is installed or not.
    printf "\nDetermining whether Vagrant needs to be installed or not.\n\n"

    if test "$(command -v vagrant)"; then
        printf "\nVagrant is installed.\n\n"
        
    else
        printf "\nVagrant is not installed.\n\n"

        # Install Virtualbox (and Homebrew, if necessary).
        printf "\nInstalling dependencies of Vagrant."
        if test "$(command -v brew)"; then
        /bin/bash -c "$(curl -fsSL \
        https://raw.githubusercontent.com/Homebrew/install/master/install.sh)" \
        && printf "\nDone installing homebrew."
        fi

        if test "$(command -v virtualbox)"; then
        brew cask install virtualbox \
        && printf "\nDone installing virtualbox."
        fi

        # Install Vagrant on MacOS.
        brew cask install vagrant && \
        brew cask install vagrant-manager && \
        printf "\nDone installing vagrant.\n"

    fi

    # Look for an existing singularity .sif file.
    if [ -f "singularity.sif" ]; then
        printf "\nSingularity .sif file found.\n\n"
        
    else
        printf "\nSingularity .sif file not found.\n\n"

        ## Proceed with building a .sif file.
        #printf "\tBuilding singularity container image file.\n\n"
        #/bin/bash singularity_build_on_mac.sh

        # Proceed with pulling a .sif file.
        printf "\nDownloading singularity container image file.\n\n"
        /bin/bash singularity_pull.sh
    fi


elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    printf "\nDetected Linux."
    # Do something under GNU/Linux platform

    if test "$(command -v singularity)"; then
        printf "\nSingularity is installed.\n\n"

        # Look for an existing singularity .sif file.
        if [ -f "singularity.sif" ]; then
            printf "\nSingularity .sif file found.\n\n"
            
        else
            printf "\nSingularity .sif file not found.\n\n"

            ## Proceed with building a .sif file.
            #printf "\tBuilding singularity container image file.\n\n"
            #/bin/bash singularity_build.sh

            # Proceed with downloading a .sif file.
            printf "\nDownloading singularity container image file using the singularity_pull.sh script.\n\n"
            /bin/bash singularity_pull.sh
        fi
        
    else
        printf "\nSingularity is not installed.\n\n"

        # Ask user if they want to install virtualbox and vagrant for running
        # a vm with singularity in it (this option will work well for desktops
        # and laptops).
        printf "\nProcedures for installing singularity on linux vary. If
        singularity is not installed on your system, consult the singularity
        documentation, install singularity, and re-run this script. This will
        not work on clusters where you do not have administrative privileges,
        so have your system administrator install singularity if necessary.

            https://sylabs.io/guides/3.0/user-guide/installation.html

        Alternatively, if you are working with linux on a personal computer
        (laptop or desktop machine) then you may find it easier to run a
        virtual machine with singularity and dependencies already installed,
        this requires VirtualBox and Vagrant, which can be installed
        automatically.
        "
        # Prompt user.
        while true; do
            read -p "Do you wish to run singularity using VirtualBox and Vagrant? (yes or no)" yn
            case $yn in
                [Yy]* ) break;;
                [Nn]* ) exit;;
                * ) echo "Please answer yes or no.";;
            esac
        done

        # Optionally install virtualbox and vagrant, similarly to how it is
        # done on MacOS.
        # This code is run if user selects yes when prompted above.

        # Check whether Vagrant is installed or not.
        printf "\nDetermining whether Vagrant needs to be installed or not.\n\n"

        if test "$(command -v vagrant)"; then
            printf "\nVagrant is installed.\n\n"
            
        else
            printf "\nVagrant is not installed.\n\n"

            # Install Virtualbox.
            printf "\nInstalling dependencies of Vagrant."
                if test "$(command -v virtualbox)"; then

                sudo apt-get install virtualbox \
                && printf "\nDone installing virtualbox."

            fi

            # Upgrade apt.
            sudo apt upgrade && \
            #sudo apt install vagrant
            # Download vagrant.
            curl -O https://releases.hashicorp.com/vagrant/2.1.0/vagrant_2.1.0_x86_64.deb && \
            # Install Vagrant on Linux.
            sudo apt install ./vagrant_2.1.0_x86_64.deb && \
            # Install plugin for setting disk size in Vagrantfile.
            #vagrant plugin install vagrant-disksize && \
            printf "\nDone installing vagrant."


        fi

        # Proceed with downloading a .sif file.
        printf "\nDownloading singularity container image file using the singularity_pull.sh script.\n\n"
        /bin/bash singularity_pull.sh

    fi

elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW32_NT" ]; then
    printf "\nDetected Windows."
    # Do something under 32 bits Windows NT platform
    echo "AMOEBAE has not been tested on Windows, please contact the developer."
    exit 1

elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW64_NT" ]; then
    printf "\nDetected Windows."
    # Do something under 64 bits Windows NT platform
    echo "AMOEBAE has not been tested on Windows, please contact the developer."

    exit 1
fi



#######################################################
# Run test code using pytest. This runs all the tests defined in files in the
# tests subdirectory.
#printf "\nTesting AMOEBAE source code using PyTest..." &&\
#docker run --rm -it \
#    -v ${PWD}:/opt/amoebae \
#    $DI /bin/bash -c "pytest"


#######################################################
# Report how much time it took for this script to run.
echo ""
ELAPSED="Setting up AMOEBAE took the following amount of time: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED


