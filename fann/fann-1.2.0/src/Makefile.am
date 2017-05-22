SUBDIRS = include

lib_LTLIBRARIES = libfloatfann.la libdoublefann.la libfixedfann.la libfann.la

AM_LDFLAGS = -version-info 3:0:2
AM_CFLAGS = -D_REENTRANT

libfloatfann_la_SOURCES = floatfann.c
libdoublefann_la_SOURCES = doublefann.c
libfixedfann_la_SOURCES = fixedfann.c
libfann_la_SOURCES = fann.c fann_io.c fann_train.c fann_train_data.c fann_options.c fann_error.c
