
# macro definitions

PROG = /usr/local/bin/get_passes
LIBS = -lm
OPTIM = 

OBJS =\
cent_gravity.o\
count_pixels.o\
do_strings.o\
fields.o\
get_app_cent.o\
get_data.o\
get_ts_data.o\
get_km_data.o\
info_content.o\
get_passes.o\
threshold.o 

# program dependencies

$(PROG):$(OBJS)
	cc $(OPTIM) -o $(PROG) $(OBJS) $(LIBS)


# inference rule

.c.o:
	cc -c $(OPTIM) $*.c
	chmod a+rw $*.o


# function dependencies

cent_gravity.o:			time_series.h
count_pixels.o:			time_series.h
do_strings.o:			time_series.h
fields.o:			time_series.h
get_app_cent.o:			time_series.h
get_data.o:			time_series.h
get_km_data.o:			km.h time_series.h
get_ts_data.o:			time_series.h
get_passes.o:			time_series.h get_passes.h
threshold.o:			time_series.h

