def options(opt):
    opt.load ('compiler_c')

def configure(conf):
    conf.define (key='PACKAGE', val='pucket')
    conf.define (key='VERSION', val='pre')

    conf.load ('compiler_c')

    conf.env.append_unique ('CFLAGS', '-std=gnu99')
    conf.env.append_unique ('CFLAGS', '-D_GNU_SOURCE')
    conf.env.append_unique ('CFLAGS', '-fopenmp')
    conf.env.append_unique ('LINKFLAGS', '-fopenmp')

    conf.check_cfg (path='xml2-config', args='--cflags --libs', package='', uselib_store='xml2')
    conf.check_cc (lib='xml2', header_name='libxml/parser.h', function_name='xmlParseFile', use='xml2')
    conf.check_cfg (package='libpng', uselib_store='png', args=['--cflags', '--libs'], msg='Checking for library png')
    conf.check_cc (lib='amdlibm', header_name='amdlibm.h', mandatory=False, define_name='HAVE_AMDLIBM', uselib_store='amdlibm')

    # does not work
    conf.check_cc (function_name='__builtin_prefetch', define_name='HAVE_BUILTIN_PREFETCH')
    conf.write_config_header ('config.h')

def build(bld):
    bld.program (features='c cprogram', source='flam3.c parser.c variations.c interpolation.c palettes.c png.c random.c rect.c main.c genome.c palettes_builtin.c', target='pucket', use='xml2 png amdlibm')

