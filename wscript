def options(opt):
    opt.load ('compiler_c')

def configure(conf):
    conf.define (key='VERSION', val='pre')
    conf.load ('compiler_c')

    conf.env.append_unique ('CFLAGS', '-std=gnu99')
    conf.env.append_unique ('CFLAGS', '-D_GNU_SOURCE')

    conf.check_cfg (path='xml2-config', args='--cflags --libs', package='', uselib_store='xml2')
    conf.check_cc (lib='xml2', header_name='libxml/parser.h', function_name='xmlParseFile', use='xml2')
    conf.check_cc (lib='pthread', uselib_store='pthread')
    conf.check_cfg (package='libpng', uselib_store='png', args=['--cflags', '--libs'], msg='Checking for library png')
    conf.check_cc (lib='amdlibm', header_name='amdlibm.h', mandatory=False, define_name='HAVE_AMDLIBM', uselib_store='amdlibm')

    # does not work
    #conf.check_cc (function_name='__builtin_ia32_rdrand64_step', define_name='HAVE_RDRAND64')
    conf.write_config_header ('config.h')

def build(bld):
    bld.program (features='c cprogram', source='flam3.c filters.c parser.c variations.c interpolation.c palettes.c png.c random.c rect.c main.c', target='vlam3', use='xml2 png amdlibm pthread')

