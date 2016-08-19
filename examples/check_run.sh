#/bin/bash

if [ $# -eq 0 ]; then
    echo "Usage: $0 directory "; exit 1
fi

BINPATH=$1

if [ -e "dns.out" ]; then
   echo -e "\033[1;31mFailed \033[0m[dns.out exists]."; exit 1
else

if [ -e "dns.ini" ]; then

#PreProcessing
    $BINPATH/inigrid.x
    if [ $? = 0 ]; then
        $BINPATH/inirand.x
        if [ $? = 0 ]; then
	    $BINPATH/iniscal.x
	    if [ $? = 0 ]; then
		$BINPATH/iniflow.x
		if [ $? = 0 ]; then
		    $BINPATH/inipart.x

#Simulation
		    if [ $? = 0 ]; then
			LIST=`ls *.ics*`; for FILE in $LIST; do mv $FILE ${FILE/ics/0}; done

			$BINPATH/dns.x
  			if [[ $? = 0 && ! -e "dns.err" ]]; then
			    diff dns.out dns.out.ref > /dev/null 2>&1
   			    if [ $? = 0 ]; then
				grep -i " nan " avg* > /dev/null 2>&1
				if [ $? = 1 ]; then
				    echo -e "\033[1;32mPassed\033[0m."
				else
				    echo -e "\033[1;31mFailed \033[0m[NaN in averages]."; exit 10
				fi
			    else
				echo -e "\033[1;31mFailed \033[0m[dns.out]."; exit 9
			    fi
			else
    			    echo -e "\033[1;31mFailed \033[0m[dns]."; exit 8
			fi
		    else
			echo -e "\033[1;31mFailed \033[0m[inipart]."; exit 7
		    fi
		else
		    echo -e "\033[1;31mFailed \033[0m[iniflow]."; exit 6
		fi
            else
		echo -e "\033[1;31mFailed \033[0m[[iniscal]."; exit 5
            fi
	else
            echo -e "\033[1;31mFailed \033[0m[[inirand]."; exit 4
	fi
    else
	echo -e "\033[1;31mFailed \033[0m[[inigrid]."; exit 3
    fi      

else
    echo -e "\033[1;31mFailed \033[0m[[dns.ini]."; exit 2
fi 

fi
