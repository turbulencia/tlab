#ifdef USE_BLOCK
     OPEN(LOC_UNIT_ID,file=name,status=LOC_STATUS,form='unformatted',block=0)
#else
#ifdef USE_ACCESS_STREAM
     OPEN(LOC_UNIT_ID,file=name,status=LOC_STATUS,form='unformatted',access='stream')
#else
#ifdef USE_BUFLESS
     OPEN(LOC_UNIT_ID,file=name,status=LOC_STATUS,form='unformatted',type='bufless')
#else
#ifdef USE_BINARY
     OPEN(LOC_UNIT_ID,file=name,status=LOC_STATUS,form='binary')
#else
     OPEN(LOC_UNIT_ID,file=name,status=LOC_STATUS,form='unformatted')
#endif
#endif
#endif
#endif
