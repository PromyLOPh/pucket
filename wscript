def options(opt):
    opt.load ('compiler_c')

def configure(conf):
    conf.define (key='VERSION', val='pre')
    conf.load ('compiler_c')

    conf.env.append_unique ('CFLAGS', '-std=gnu99')
    conf.env.append_unique ('CFLAGS', '-mrdrnd')

    conf.check_cfg (path='xml2-config', args='--cflags --libs', package='', uselib_store='xml2')
    conf.check_cc (lib='xml2', header_name='libxml/parser.h', function_name='xmlParseFile', use='xml2')
    conf.check_cc (lib='pthread', uselib_store='pthread')
    conf.check_cc (lib='jpeg', uselib_store='jpeg')
    conf.check_cfg (package='libpng', uselib_store='png', args=['--cflags', '--libs'], msg='Checking for library png')
    conf.check_cc (lib='amdlibm', header_name='amdlibm.h', mandatory=False, define_name='HAVE_AMDLIBM', uselib_store='amdlibm')

    # does not work
    #conf.check_cc (function_name='__builtin_ia32_rdrand64_step', define_name='HAVE_RDRAND64')
    conf.write_config_header ('config.h')

def build(bld):
    bld.stlib (features='c cstlib', source='flam3.c filters.c parser.c variations.c interpolation.c palettes.c jpeg.c png.c xorshift.c docstring.c', target='libflam3', use='xml2 png jpeg pthread', includes='.')
    bld.program (features='c cprogram', source='flam3-render.c', target='flam3-render', use='libflam3 xml2 jpeg png amdlibm pthread', includes='.')
    bld.program (features='c cprogram', source='flam3-genome.c', target='flam3-genome', use='libflam3 xml2 png amdlibm pthread', includes='.')
    bld.program (features='c cprogram', source='flam3-animate.c', target='flam3-animate', use='libflam3 xml2 png amdlibm pthread', includes='.')
    #bld.program (features='c cprogram', source='flam3-convert.c', target='flam3-convert', use='libflam3 xml2 png amdlibm pthread', includes='.')

