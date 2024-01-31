#!/bin/sh
######################################################################
# Prepare a Job to be submited to the system background queues. 
# Queue manager depends on site and system.
######################################################################
export TOOLS_HOME=$HOME/bin
export CONFIG_FILE=$TOOLS_HOME/CONFIG.default

QSEND_TEST="no"

. $CONFIG_FILE

# Output Configuration
echo ""
echo "Qsend version 2.0"
echo ""
echo "Site    = $site"
echo "Account = $account"
echo "Queue   = $queue"
echo "HPC mode= $hpc"
echo ""
echo "Current Directory = `pwd`"
echo ""

# Timestamp of this submission
TIMESTAMP=`date +%Y%m%d-%H%M%S`

# Temporary file
QSEND_TMP_SCRIPT="dns.nqs.sh-$TIMESTAMP"

# Default values
QSEND_NAME="dns-run"
QSEND_JOBDIR=`pwd`
QSEND_JOBDIR=${QSEND_JOBDIR#"$WORK"/}
if [ $hpc = "parallel" ]; then
   QSEND_NP=32
else
   QSEND_NP=1
fi
QSEND_MEM=512
QSEND_TIME=14400
QSEND_MAXITER=2000
QSEND_STEP=500
QSEND_SCRIPT="dns.nqs.sh.$site.$system"
QSEND_QUEUE="$queue"
if [ $system = "jugene" ]; then
   QSEND_CPUMODE="smp"
fi
if [ $system = "juqueen" ]; then
   QSEND_COMMMODE=blocking
   QSEND_RANKSPERNODE=16
fi
if [ $system = "juwels"  ]; then
   QSEND_COMMMODE=blocking
   QSEND_RANKSPERNODE=48
   QSEND_THREADSPERRANK=1
   QSEND_QUEUE=batch
fi
if [ $system = "jureca" ]; then 
   QSEND_COMMODE=blocking 
   QSEND_RANKSPERNODE=24
   QSEND_THREADSPERRANK=1
fi
if [ $system = "supermuc" ]; then
   QSEND_COMMMODE=blocking
   QSEND_RANKSPERNODE=48
   QSEND_THREADSPERRANK=1
fi

RUNMODE="simulation"

# Check if options passed
if [ $# -eq 0 ]; then

###################################################
# Getting interactive information
###################################################
   interactive="yes"

   echo "Enter Job Name [$QSEND_NAME] :"
   read value
   if [ -n "$value" ]; then
      QSEND_NAME=$value
   fi

   echo "Enter Job Directory [$QSEND_JOBDIR] :"
   read value
   if [ -n "$value" ]; then
      QSEND_JOBDIR=$value
   fi
   # read QSEND_JOBDIR
   # if [ -z "$QSEND_JOBDIR" ]; then
   #    echo "Error: Empty JOBDIR"
   #    exit 1
   # fi

   if [ $hpc = "parallel" ]; then
      echo "Enter Number of Tasks :"
      read value
      if [ -n "$value" ]; then
	  QSEND_NP=$value
      fi
   fi

   if [ $system = "blizzard" ]; then
       echo "Enter Memory Requested :Mega(words/bytes)"
       read value
       if [ -n "$value" ]; then
	   QSEND_MEM=$value
       fi
   fi
  
   if [ $system = "jugene" ]; then
      echo "Enter cpu mode :(vn/dual/smp)"
      read value
      if [ -n "$value" ]; then
         QSEND_CPUMODE=$value
      fi
   fi
   
   if [ $system = "juqueen" ]; then
      echo "Enter communication mode :(blocking/nonblocking)"
      read value
      if [ -n "$value" ]; then
         QSEND_COMMMODE=$value
      fi

      echo "Enter ranks per node :(16/32)"
      read value
      if [ -n "$value" ]; then
         QSEND_RANKSPERNODE=$value
      fi

   fi

   if [ $system = "juwels"  ]; then 
      echo "Enter ranks per node [$QSEND_RANKSPERNODE] :(1-48)"
      read value
      if [ -n "$value" ]; then
         QSEND_RANKSPERNODE=$value
      fi

      QSEND_THREADSPERRANK=$((48/$QSEND_RANKSPERNODE))
      echo "Enter threads per rank [$QSEND_THREADSPERRANK] :(1-48)"
      read value
      if [ -n "$value" ]; then
         QSEND_THREADSPERRANK=$value
      fi

      echo "Enter system partition [$QSEND_QUEUE] : (batch,mem192,devell)"
      read value
      if [ -n "$value" ]; then  
	  QSEND_QUEUE=$value
      fi 

   fi

   if [ $system = "supermuc" ]; then 
      echo "Enter ranks per node [$QSEND_RANKSPERNODE] :(1-48)"
      read value
      if [ -n "$value" ]; then
         QSEND_RANKSPERNODE=$value
      fi

      QSEND_THREADSPERRANK=$((48/$QSEND_RANKSPERNODE))
      echo "Enter threads per rank [$QSEND_THREADSPERRANK] :(1-48)"
      read value
      if [ -n "$value" ]; then
         QSEND_THREADSPERRANK=$value
      fi

      echo "Enter system partition [$QSEND_QUEUE] : (general,micro,test)"
      read value
      if [ -n "$value" ]; then  
	  QSEND_QUEUE=$value
      fi 

   fi

   if [ $system = "jureca" ] ; then 
      echo "Enter ranks per node [$QSEND_RANKSPERNODE] :(1-24)"
      read value
      if [ -n "$value" ]; then
         QSEND_RANKSPERNODE=$value
      fi

      QSEND_THREADSPERRANK=$((24/$QSEND_RANKSPERNODE))
      echo "Enter threads per rank [$QSEND_THREADSPERRANK] :(1-24)"
      read value
      if [ -n "$value" ]; then
         QSEND_THREADSPERRANK=$value
      fi
   fi

   echo "Enter Time Requested :(seconds)"
   read value
   if [ -n "$value" ]; then
      QSEND_TIME=$value
   fi

   echo "Enter Run Mode [$RUNMODE] :(preprocess/simulation/postprocess)"
   read value
   if [ -n "$value" ]; then
      RUNMODE=$value
   fi

   if [ $RUNMODE = "simulation" ]; then

      echo "Maximum iterations :(1000/2500/5000)"
      read value
      if [ -n "$value" ]; then
         QSEND_MAXITER=$value
      fi

      echo "Enter Rerun Step :(500)"
      read value
      if [ -n "$value" ]; then
         QSEND_STEP=$value
      fi

   fi

   echo "Enter Script to Execute ($QSEND_SCRIPT)"
   read value
   if [ -n "$value" ]; then
      QSEND_SCRIPT=$value
   fi

   echo ""

   echo "Job Name             : $QSEND_NAME"
   echo "Job Directory        : $QSEND_JOBDIR"
   if [ $hpc = "parallel" ]; then
      echo "Processors Requested : $QSEND_NP"
   fi 
   echo "Memory Requested     : $QSEND_MEM"
   echo "Time requested       : $QSEND_TIME"
   echo "Run mode             : $RUNMODE"
   if [ $RUNMODE = "simulation" ]; then
      echo "Maximum Iteration    : $QSEND_MAXITER"
      echo "Rerun Step           : $QSEND_STEP"
   fi
   if [ $system = "jugene" ]; then
      echo "CPU mode             : $QSEND_CPUMODE"
   fi
   if [ $system = "juqueen" ]; then
      echo "Communication mode   : $QSEND_COMMMODE"
      echo "Ranks per node       : $QSEND_RANKSPERNODE"
   fi
   if [ $system = "juwels" ] || [ $system = "jureca" ]  || [ $system = "supermuc" ]; then 
      echo "Partition            : $QSEND_QUEUE"
      echo "Ranks per node       : $QSEND_RANKSPERNODE"
      echo "Threads per rank     : $QSEND_THREADSPERRANK"
   fi
   echo "Script               : $QSEND_SCRIPT"

   echo ""

   echo "Is this correct (yes/no)"
   read answer
   if [ -z "$answer" ]; then
      echo "You must answer yes or no"
      exit 1
   fi

else

###################################################
# Getting non-interactive information
###################################################
  interactive="no"

  for opt in $*; do
    case $opt in
	-name)
	    shift
	    QSEND_NAME=$1
	    shift
	    ;;
	-queue)
	    shift
	    QSEND_QUEUE=$1
	    shift
	    ;; 
	-time)
	    shift
	    QSEND_TIME=$1
	    shift
	    ;;
        -maxiter) 
            shift
            QSEND_MAXITER=$1
            shift
            ;;
	-mem)
	    shift
	    QSEND_MEM=$1
	    shift
	    ;;
	-jobdir)
	    shift
	    QSEND_JOBDIR=$1
	    shift
	    ;;
	-script)
	    shift
	    QSEND_SCRIPT=$1
	    shift
	    ;;	
	-test)
	    shift
	    QSEND_TEST="yes"
	    ;;
        -step)
            shift
            QSEND_STEP=$1
            shift
            ;;
	-np)
	    shift
            QSEND_NP=$1
            shift
	    ;;
	-cpumode)
	    shift
            QSEND_CPUMODE=$1
            shift
	    ;;
	-rankspernode)
	    shift
            QSEND_RANKSPERNODE=$1
            shift
	    ;;
	-threadsperrank)
	    shift
            QSEND_THREADSPERRANK=$1
            shift
	    ;;
	-commmode)
	    shift
            QSEND_COMMMODE=$1
            shift
	    ;;
	esac
    done	
fi

if [ -z "$QSEND_JOBDIR" ]; then
   echo "Error: Empty JOBDIR"
   exit 1
fi

# Check quota 
#BALANCE=`$TOOLS_HOME/dns.account $account`

if [ -z "$BALANCE" ]; then
   BALANCE="999999999"
fi

if [ $BALANCE -lt $QSEND_TIME ]; then
   echo "You don't have enough time available."
   exit 1
fi

if [ $interactive = "no" ]; then
   answer="yes"
fi

# Creating script file and submitting
if [ $answer = "yes" ]; then
	
###################################################
# Writting information in script file
###################################################
# CRAY T3E
    if [ $system = "t3e" ]; then
	LINE='BEGIN {} 
	    { 
	    if ($2 == "-r" )
		{
		    print "#QSUB -r \"'$QSEND_NAME'\"" 
		}
	    else if ($4 == "DIRECTIVES" ) 
		{
		    print "#QSUB -l mpp_p='$QSEND_NP'"
		    print "#QSUB -l mpp_t='$QSEND_TIME'"
		    print $0
		}
	    else if ($4 == "MODIFY")
		{
		    print "JOBDIR='$QSEND_JOBDIR'"
		    print "MEM='$QSEND_MEM'"
		    print "TIME='$QSEND_TIME'"
		    print "NAME=\"'$QSEND_NAME'\""
                    print "MAXITER=\"'$QSEND_MAXITER'\""
		    print "SCRIPT='$QSEND_SCRIPT'"
                    print "STEP='$QSEND_STEP'"
                    print "NP='$QSEND_NP'"
                    print "ACCOUNT='$account'"
                    print "TIMESTAMP='$TIMESTAMP'"
		    print $0
		}
	    else
		{
		    print $0
		}
	}' 
# LINUX 
    elif [ $system = "x86" ]; then
	if [ $site = "dkrz" ]; then
	LINE='BEGIN {} 
	    { 
             if ($4 == "DIRECTIVES" ) 
		{
		    print "#SBATCH --job-name='$QSEND_NAME'"
                    print "#SBATCH --output='$QSEND_NAME'.o%j"
                    print "#SBATCH --error='$QSEND_NAME'.e%j"
                    print "#SBATCH --ntasks='$QSEND_NP'.e%j"
 
		    total_time = '$QSEND_TIME';
		    hours = int(total_time/3600);
		    total_time = total_time-hours*3600;
		    min   = int(total_time/60);
		    sec   = int(total_time-min*60);
		    printf "#SBATCH --time=%.2d:%.2d:%.2d\n", hours, min, sec
		}
	    else if ($4 == "MODIFY")
		{
		    print "JOBDIR='$QSEND_JOBDIR'"
		    print "MEM='$QSEND_MEM'"
		    print "TIME='$QSEND_TIME'"
		    print "NAME=\"'$QSEND_NAME'\""
                    print "MAXITER=\"'$QSEND_MAXITER'\""
		    print "SCRIPT='$QSEND_SCRIPT'"
                    print "STEP='$QSEND_STEP'"
                    print "NP='$QSEND_NP'"
                    print "ACCOUNT='$account'"
                    print "TIMESTAMP='$TIMESTAMP'"
		}
	    print $0
	}' 
	elif [ $site = "tacc" ]; then
	LINE='BEGIN {} 
	    { 
	    if ($2 == "-o" )
		{
		    print "#BSUB -J '$QSEND_NAME'"
		}
	    else if ($4 == "DIRECTIVES" ) 
		{
		    total_time = '$QSEND_TIME';
		    hours = int(total_time/3600);
		    total_time = total_time-hours*3600;
		    min   = int(total_time/60);
		    sec   = int(total_time-min*60);
		    
		    print "#BSUB -n " '$QSEND_NP'
		    print "#BSUB -W " hours ":" min 
		}
	    else if ($4 == "MODIFY")
		{
		    print "JOBDIR='$QSEND_JOBDIR'"
		    print "MEM='$QSEND_MEM'"
		    print "TIME='$QSEND_TIME'"
		    print "NAME=\"'$QSEND_NAME'\""
                    print "MAXITER=\"'$QSEND_MAXITER'\""
		    print "SCRIPT='$QSEND_SCRIPT'"
                    print "STEP='$QSEND_STEP'"
                    print "NP='$QSEND_NP'"
                    print "ACCOUNT='$account'"
                    print "TIMESTAMP='$TIMESTAMP'"
		}
	    print $0
	}' 
        elif [ $site = "rwth" ]; then
	LINE='BEGIN {} 
	    { 
	    if ($4 == "DIRECTIVES" ) 
		{
		    total_time = '$QSEND_TIME';
		    hours = int(total_time/3600);
		    total_time = total_time-hours*3600;
		    min   = int(total_time/60);
		    sec   = int(total_time-min*60);

		    print "#$ -N '$QSEND_NAME'"
		    print "#$ -o qsub_out.log-'$TIMESTAMP'"
		    print "#$ -e qsub_err.log-'$TIMESTAMP'"
		    print "#$ -l h_vmem='$QSEND_MEM'M"
		    print "#$ -l h_rt=" hours ":" min ":" sec
                    print "#$ -pe mpi* '$QSEND_NP'"
		}
	    else if ($4 == "MODIFY")
		{
		    print "JOBDIR='$QSEND_JOBDIR'"
		    print "MEM='$QSEND_MEM'"
		    print "TIME='$QSEND_TIME'"
		    print "NAME=\"'$QSEND_NAME'\""
                    print "MAXITER=\"'$QSEND_MAXITER'\""
		    print "SCRIPT='$QSEND_SCRIPT'"
                    print "STEP='$QSEND_STEP'"
                    print "NP='$QSEND_NP'"
                    print "ACCOUNT='$account'"
                    print "TIMESTAMP='$TIMESTAMP'"
                    print "RUNMODE='$RUNMODE'"
		}
	    print $0
	}' 
        else
	LINE='BEGIN {} 
	    { 
	    if ($4 == "DIRECTIVES" ) 
		{
		    print "#$ -N '$QSEND_NAME'"
		    print "#$ -o qsub_out.log-'$TIMESTAMP'"
		    print "#$ -e qsub_err.log-'$TIMESTAMP'"
		    print "#$ -l h_vmem='$QSEND_MEM'M"
		}
	    else if ($4 == "MODIFY")
		{
		    print "JOBDIR='$QSEND_JOBDIR'"
		    print "MEM='$QSEND_MEM'"
		    print "TIME='$QSEND_TIME'"
		    print "NAME=\"'$QSEND_NAME'\""
                    print "MAXITER=\"'$QSEND_MAXITER'\""
		    print "SCRIPT='$QSEND_SCRIPT'"
                    print "STEP='$QSEND_STEP'"
                    print "NP='$QSEND_NP'"
                    print "ACCOUNT='$account'"
                    print "TIMESTAMP='$TIMESTAMP'"
                    print "RUNMODE='$RUNMODE'"
		}
	    print $0
	}' 
        fi
# x86 mistral
    elif [ $system = "mistral" ]; then
	LINE='BEGIN {} 
	    { 
             if ($4 == "DIRECTIVES" ) 
		{
		    print "#SBATCH --job-name='$QSEND_NAME'"
                    print "#SBATCH --output=ll_out.%j-'$TIMESTAMP'"
                    print "#SBATCH --error=ll_err.%j-'$TIMESTAMP'"
                    print "#SBATCH --ntasks='$QSEND_NP'"
 
		    total_time = '$QSEND_TIME';
		    hours = int(total_time/3600);
		    total_time = total_time-hours*3600;
		    min   = int(total_time/60);
		    sec   = int(total_time-min*60);
		    printf "#SBATCH --time=%.2d:%.2d:%.2d\n", hours, min, sec
		}
	    else if ($4 == "MODIFY")
		{
		    print "JOBDIR='$QSEND_JOBDIR'"
		    print "MEM='$QSEND_MEM'"
		    print "TIME='$QSEND_TIME'"
		    print "NAME=\"'$QSEND_NAME'\""
                    print "MAXITER=\"'$QSEND_MAXITER'\""
		    print "SCRIPT='$QSEND_SCRIPT'"
                    print "STEP='$QSEND_STEP'"
                    print "NP='$QSEND_NP'"
                    print "ACCOUNT='$account'"
                    print "TIMESTAMP='$TIMESTAMP'"
                    print "RUNMODE='$RUNMODE'"
		}
	    print $0
	}' 
    elif [ $system = "juwels" ] || [ $system = "jureca" ]  || [ $system = "supermuc" ]; then
	LINE='BEGIN {} 
	    { 
             if ($4 == "DIRECTIVES" ) 
		{
		    print "#SBATCH --job-name='$QSEND_NAME'"
                    print "#SBATCH --output=slurm_out.%j-'$TIMESTAMP'"
                    print "#SBATCH --error=slurm_err.%j-'$TIMESTAMP'"

		    total_pe = '$QSEND_NP';
		    rpn      = '$QSEND_RANKSPERNODE';
                    nodes = int((total_pe-1)/rpn) + 1;

                    print "#SBATCH --partition='$QSEND_QUEUE'" 
                    print "#SBATCH --nodes=" nodes
                    print "#SBATCH --ntasks-per-node='$QSEND_RANKSPERNODE'"
                    if ('$QSEND_THREADSPERRANK' >= 2)
                       print "#SBATCH --cpus-per-task='$QSEND_THREADSPERRANK'"
 
		    total_time = '$QSEND_TIME';
		    hours = int(total_time/3600);
		    total_time = total_time-hours*3600;
		    min   = int(total_time/60);
		    sec   = int(total_time-min*60);
		    printf "#SBATCH --time=%.2d:%.2d:%.2d\n", hours, min, sec

		}
	    else if ($4 == "MODIFY")
		{
		    print "JOBDIR='$QSEND_JOBDIR'"
                    print "QUEUE='$QSEND_QUEUE'"
		    print "MEM='$QSEND_MEM'"
		    print "TIME='$QSEND_TIME'"
		    print "NAME=\"'$QSEND_NAME'\""
                    print "MAXITER=\"'$QSEND_MAXITER'\""
		    print "SCRIPT='$QSEND_SCRIPT'"
                    print "STEP='$QSEND_STEP'"
                    print "NP='$QSEND_NP'"
                    print "ACCOUNT='$account'"
                    print "TIMESTAMP='$TIMESTAMP'"
                    print "RUNMODE='$RUNMODE'"
                    print "RANKSPERNODE='$QSEND_RANKSPERNODE'"
                    print "THREADSPERRANK='$QSEND_THREADSPERRANK'"
                    print "COMMMODE='$QSEND_COMMMODE'"
		}
	    print $0
	}' 
# SUN
    elif [ $system = "sun" ]; then
	LINE='BEGIN {} 
	    { 
	    if ($4 == "DIRECTIVES" ) 
		{
		    total_time = '$QSEND_TIME';
		    hours = int(total_time/3600);
		    total_time = total_time-hours*3600;
		    min   = int(total_time/60);
		    sec   = int(total_time-min*60);

		    print "#$ -N '$QSEND_NAME'"
		    print "#$ -o qsub_out.log-'$TIMESTAMP'"
		    print "#$ -e qsub_err.log-'$TIMESTAMP'"
		    print "#$ -l h_vmem='$QSEND_MEM'M"
		    print "#$ -l h_rt=" hours ":" min ":" sec
                    print "#$ -pe mpi_sunos* '$QSEND_NP'"
		}
	    else if ($4 == "MODIFY")
		{
		    print "JOBDIR='$QSEND_JOBDIR'"
		    print "MEM='$QSEND_MEM'"
		    print "TIME='$QSEND_TIME'"
		    print "NAME=\"'$QSEND_NAME'\""
                    print "MAXITER=\"'$QSEND_MAXITER'\""
		    print "SCRIPT='$QSEND_SCRIPT'"
                    print "STEP='$QSEND_STEP'"
                    print "NP='$QSEND_NP'"
                    print "ACCOUNT='$account'"
                    print "TIMESTAMP='$TIMESTAMP'"
                    print "RUNMODE='$RUNMODE'"
		}
	    print $0
	}' 
# HITACHI SR8000
    elif [ $system = "sr8000" ]; then
	LINE='BEGIN {} 
	    { 
	    if ($2 == "-r" )
		{
		    print "#@$ -r '$QSEND_NAME'"
		}
	    else if ($4 == "DIRECTIVES" ) 
		{
		    print "#@$-r '$QSEND_NAME'"

		    total_time = '$QSEND_TIME';
		    hours = int(total_time/3600);
		    total_time = total_time-hours*3600;
		    min   = int(total_time/60);
		    sec   = int(total_time-min*60);
		    
		    print "#@$-lE " hours ":" min ":" sec

		    total_pe = '$QSEND_NP'
		    nodes    = int((total_pe+7)/8)

		    print "#@$-N " nodes
		}
	    else if ($4 == "MODIFY")
		{
		    print "JOBDIR='$QSEND_JOBDIR'"
		    print "MEM='$QSEND_MEM'"
		    print "TIME='$QSEND_TIME'"
		    print "NAME=\"'$QSEND_NAME'\""
                    print "MAXITER=\"'$QSEND_MAXITER'\""
		    print "SCRIPT='$QSEND_SCRIPT'"
                    print "STEP='$QSEND_STEP'"
                    print "NP='$QSEND_NP'"
                    print "ACCOUNT='$account'"
                    print "TIMESTAMP='$TIMESTAMP'"

		    total_pe = '$QSEND_NP'
		    nodes    = int((total_pe+7)/8)
		    
		    print "NN=" nodes
		}
	    print $0
	}' 
# IBM Power6 LoadLeveler
    elif [[ $system = "jump" || $system = "blizzard" ]]; then
	LINE='BEGIN {} 
	    { 
	    if ($4 == "DIRECTIVES" ) 
		{
		    total_time = '$QSEND_TIME';
		    hours = int(total_time/3600);
		    total_time = total_time-hours*3600;
		    min   = int(total_time/60);
		    sec   = int(total_time-min*60);
		    
		    total_pe = '$QSEND_NP'
		    nodes    = int((total_pe-1)/32) + 1

		    print "#@ job_name = '$QSEND_NAME'"
		    print "#@ output = ll_out.$(jobid).log-'$TIMESTAMP'"
		    print "#@ error= ll_err.$(jobid).log-'$TIMESTAMP'"
		    print "#@ node = " nodes
                    print "#@ resources = ConsumableMemory('$QSEND_MEM' mb)"
                    print "#@ tasks_per_node = 32"
		    print "#@ wall_clock_limit = " hours ":" min ":" sec
		}
	    else if ($4 == "MODIFY")
		{
		    print "JOBDIR='$QSEND_JOBDIR'"
		    print "MEM='$QSEND_MEM'"
		    print "TIME='$QSEND_TIME'"
		    print "NAME=\"'$QSEND_NAME'\""
                    print "MAXITER=\"'$QSEND_MAXITER'\""
		    print "SCRIPT='$QSEND_SCRIPT'"
                    print "STEP='$QSEND_STEP'"
                    print "NP='$QSEND_NP'"
                    print "ACCOUNT='$account'"
                    print "TIMESTAMP='$TIMESTAMP'"
                    print "RUNMODE='$RUNMODE'"
		}
	    print $0
	}' 
# IBM JUGENE
    elif [ $system = "jugene" ]; then
	LINE='BEGIN {} 
	    { 
	    if ($4 == "DIRECTIVES" ) 
		{
		    total_time = '$QSEND_TIME';
		    hours = int(total_time/3600);
		    total_time = total_time-hours*3600;
		    min   = int(total_time/60);
		    sec   = int(total_time-min*60);
		    
		    total_pe = '$QSEND_NP';
                    if      ( "'$QSEND_CPUMODE'" == "vn" ) 
                       { nodes = int((total_pe-1)/4) + 1 }
                    else if ( "'$QSEND_CPUMODE'" == "dual" )
                       { nodes = int((total_pe-1)/2) + 1 }
                    else if ( "'$QSEND_CPUMODE'" == "smp" )
		       { nodes = total_pe }
                    fi

                    print  "#@ job_name = '$QSEND_NAME'"
		    print  "#@ output = ll_out.$(jobid).log-'$TIMESTAMP'"
		    print  "#@ error= ll_err.$(jobid).log-'$TIMESTAMP'"
 		    print  "#@ bg_size = " nodes
		    printf "#@ wall_clock_limit = %.2d:%.2d:%.2d\n", hours, min, sec
		}
	    else if ($4 == "MODIFY")
		{
		    print "JOBDIR='$QSEND_JOBDIR'"
		    print "MEM='$QSEND_MEM'"
		    print "TIME='$QSEND_TIME'"
		    print "NAME=\"'$QSEND_NAME'\""
                    print "MAXITER=\"'$QSEND_MAXITER'\""
		    print "SCRIPT='$QSEND_SCRIPT'"
                    print "STEP='$QSEND_STEP'"
                    print "NP='$QSEND_NP'"
                    print "ACCOUNT='$account'"
                    print "TIMESTAMP='$TIMESTAMP'"
                    print "RUNMODE='$RUNMODE'"
                    print "CPUMODE='$QSEND_CPUMODE'"
		}
	    print $0
	}' 
# IBM JUQUEEN
    elif [ $system = "juqueen" ]; then
	LINE='BEGIN {} 
	    { 
	    if ($4 == "DIRECTIVES" ) 
		{
		    total_time = '$QSEND_TIME';
		    hours = int(total_time/3600);
		    total_time = total_time-hours*3600;
		    min   = int(total_time/60);
		    sec   = int(total_time-min*60);
		    
		    total_pe = '$QSEND_NP';
		    rpn      = '$QSEND_RANKSPERNODE';
                    nodes = int((total_pe-1)/rpn) + 1;

                    print  "#@ job_name = '$QSEND_NAME'"
		    print  "#@ output = ll_out.$(jobid).log-'$TIMESTAMP'"
		    print  "#@ error= ll_err.$(jobid).log-'$TIMESTAMP'"
 		    print  "#@ bg_size = " nodes
		    printf "#@ wall_clock_limit = %.2d:%.2d:%.2d\n", hours, min, sec
		}
	    else if ($4 == "MODIFY")
		{
		    print "JOBDIR='$QSEND_JOBDIR'"
		    print "MEM='$QSEND_MEM'"
		    print "TIME='$QSEND_TIME'"
		    print "NAME=\"'$QSEND_NAME'\""
                    print "MAXITER=\"'$QSEND_MAXITER'\""
		    print "SCRIPT='$QSEND_SCRIPT'"
                    print "STEP='$QSEND_STEP'"
                    print "NP='$QSEND_NP'"
                    print "ACCOUNT='$account'"
                    print "TIMESTAMP='$TIMESTAMP'"
                    print "RUNMODE='$RUNMODE'"
                    print "RANKSPERNODE='$QSEND_RANKSPERNODE'"
                    print "COMMMODE='$QSEND_COMMMODE'"
		}
	    print $0
	}' 
# CRAY T90 Mostly
    else
	LINE='BEGIN {i=0} 
	    { 
	    if ($2 == "-lM")
		{
		    print "#QSUB -lM '$QSEND_MEM'MW" 
		    i = 1
		}
	    else if ($2 == "-lT" )
		{
		    print "#QSUB -lT '$QSEND_TIME'" 
		}
	    else if ($2 == "-r" )
		{
		    print "#QSUB -r \"'$QSEND_NAME'\"" 
		}
	    else if ($4 == "DIRECTIVES" ) 
		{
		    if ( i == 0 )
			{
			print "#QSUB -lM '$QSEND_MEM'MW"
			}
		    print $0
		}
	    else if ($4 == "MODIFY")
		{
		    print "JOBDIR='$QSEND_JOBDIR'"
		    print "MEM='$QSEND_MEM'"
		    print "TIME='$QSEND_TIME'"
		    print "NAME=\"'$QSEND_NAME'\""
                    print "MAXITER=\"'$QSEND_MAXITER'\""
		    print "SCRIPT='$QSEND_SCRIPT'"
                    print "STEP='$QSEND_STEP'"
		    print "NP=1"
                    print "ACCOUNT='$account'"
                    print "TIMESTAMP='$TIMESTAMP'"
		    print $0
		}
	    else
		{
		    print $0
		}
	}'
    fi

# Write info
    awk "$LINE" $TOOLS_HOME/$QSEND_SCRIPT > $QSEND_TMP_SCRIPT

    if [ $interactive = "yes" ]; then

	echo ""

	echo "Would you like to review the script (yes/no)"
	read answer
	
	if [ $answer != "no" ]; then
	    vi $QSEND_TMP_SCRIPT
	fi

	echo ""

	echo "Would you like to submit the job (yes/no)"
	read answer

    fi

    if [ $interactive = "no" ]; then
       answer="yes"
    fi

###################################################
# Sumitting script file
###################################################
    if [ $QSEND_TEST = "no" ]; then
	if [ $answer = "yes" ]; then
	    echo "Submiting Job $QSEND_NAME"
	    chmod +x $QSEND_TMP_SCRIPT
	    if [ $site = "tacc" ]; then
	       bsub < $QSEND_TMP_SCRIPT
            elif [ $site = "rwth" ]; then
	       qsub $queue $QSEND_TMP_SCRIPT
	    else
	       sbatch $QSEND_TMP_SCRIPT
	    fi
	fi
    fi
fi

