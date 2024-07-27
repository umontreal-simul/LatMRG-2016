#!/usr/bin/env python
# coding: utf-8

top = '.'
out = 'build'

from waflib import Utils

import imp
def waftool(name):
    return imp.load_module('waf_' + name, *imp.find_module(name, ['./latcommon/waftools']))

version = waftool('version')
compiler = waftool('compiler')
deps = waftool('deps')

def options(ctx):
    ctx.recurse('latcommon')
    ctx.add_option('--build-docs', action='store_true', default=False, help='build documentation')
    ctx.add_option('--testu01', action='store', help='prefix under which TestU01 is installed')
    ctx.add_option('--gmp', action='store', help='prefix under which GMP is installed')

def configure(ctx):
    ctx.recurse('latcommon')
    if not ctx.env.LATCOMMON_SUFFIX:
        ctx.fatal('Please specify the NTL number types with --ntltypes (try ./waf --help)')

    if ctx.options.testu01:
        deps.add_deps_path(ctx, 'TESTU01', ctx.options.testu01)
    if ctx.options.gmp:
        deps.add_deps_path(ctx, 'GMP', ctx.options.gmp)

    ctx.version_file()

    ctx_check = deps.shared_or_static(ctx, ctx.check)

    # GMP
    ctx_check(features='cxx cxxprogram', header_name='gmp.h')
    ctx_check(features='cxx cxxprogram', lib='gmp', uselib_store='GMP')

    # TestU01
    ctx_check(features='cxx cxxprogram', header_name='ulcg.h')
    ctx_check(features='cxx cxxprogram', lib='testu01', uselib_store='TESTU01', use='GMP')

    # mylib (part of TestU01)
    ctx_check(features='cxx cxxprogram', header_name='num.h')
    ctx_check(features='cxx cxxprogram', lib='mylib', uselib_store='MYLIB')

    # Doxygen
    if ctx.options.build_docs:
        ctx.env.BUILD_DOCS = True
        if not ctx.find_program('doxygen', var='DOXYGEN', mandatory=False):
            ctx.fatal('Doxygen is required for building documentation.\n' +
                      'Get it from http://www.stack.nl/~dimitri/doxygen/')

    # version
    ctx.define('LATMRG_VERSION', ctx.set_version())
    ctx.msg("Setting LatMRG version", version.VERSION)

    # build variants
    env = ctx.env.derive()
    env.detach()

    # release (default)
    ctx.env.append_unique('CXXFLAGS', ['-O2'])
    ctx.define('NDEBUG', 1)

    ctx.setenv('debug', env)
    ctx.env.append_unique('CXXFLAGS', ['-g'])
    ctx.define('DEBUG', 1)


def distclean(ctx):
    ctx.recurse('latcommon')
    verfile = ctx.path.find_node('VERSION')
    if verfile:
        verfile.delete()
    from waflib import Scripting
    Scripting.distclean(ctx)


def build(ctx):

    if ctx.variant:
        print("Building variant `%s'" % (ctx.variant,))

    ctx.recurse('latcommon')
    ctx.recurse('src')
    ctx.recurse('progs')
    ctx.recurse('test')
    if ctx.env.BUILD_DOCS:
        ctx.recurse('doc')


# build variants

from waflib.Build import BuildContext
class debug(BuildContext):
    cmd = 'debug'
    variant = 'debug'
