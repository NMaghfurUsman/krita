/* Generic functions for reading XCF files
 *
 * This file was written by Henning Makholm <henning@makholm.net>
 * It is hereby in the public domain.
 *
 * In jurisdictions that do not recognise grants of copyright to the
 * public domain: I, the author and (presumably, in those jurisdictions)
 * copyright holder, hereby permit anyone to distribute and use this code,
 * in source code or binary form, with or without modifications. This
 * permission is world-wide and irrevocable.
 *
 * Of course, I will not be liable for any errors or shortcomings in the
 * code, since I give it away without asking any compenstations.
 *
 * If you use or distribute this code, I would appreciate receiving
 * credit for writing it, in whichever way you find proper and customary.
 */

#include "xcftools.h"
#include <string.h>
#include <errno.h>
#include <limits.h>
#ifdef HAVE_ICONV
# include <iconv.h>
#elif !defined(ICONV_CONST)
# define ICONV_CONST const
#endif

uint8_t *xcf_file = 0 ;
size_t xcf_length ;
int use_utf8 = 0 ;

int
xcfOffset(uint32_t addr,int spaceafter, uint32_t* apparent)
{
  if (!apparent) {
      return XCF_ERROR;
  }
  if(xcfCheckspace(addr,4,"(xcfOffset)") != XCF_OK) {
      return XCF_ERROR;
  }
  *apparent = xcfL(addr);
  if (xcfCheckspace(*apparent,spaceafter,
                "Too large offset (%" PRIX32 ") at position %" PRIX32,
                                *apparent,addr) != XCF_OK) {
      return XCF_ERROR;
  }
  return XCF_OK;
}

int
xcfNextprop(uint32_t *master,uint32_t *body, PropType *typeOut)
{
  int response;

  if (typeOut ==  0) {
      return XCF_ERROR;
  }

  uint32_t ptr, length, total, minlength ;
  PropType type ;
  ptr = *master ;
  if ((response = xcfCheckspace(ptr,8,"(property header)")) != XCF_OK) {
      return XCF_ERROR;
  }
  type = xcfL(ptr);
  length = xcfL(ptr+4);
  *body = ptr+8 ;

  switch(type) {
  case PROP_COLORMAP:
    {
      uint32_t ncolors ;
      if ((response = xcfCheckspace(ptr+8,4,"(colormap length)")) != XCF_OK) {
          return XCF_ERROR;
      }
      ncolors = xcfL(ptr+8) ;
      if( ncolors > 256 ) {
        FatalBadXCF("Colormap has %" PRIu32 " entries",ncolors);
        return XCF_ERROR;
      }

      /* Surprise! Some older version of the Gimp computed the wrong length
       * word, and the _reader_ always just reads three bytes per color
       * and ignores the length tag! Duplicate this so we too can read
       * the buggy XCF files.
       */
      length = minlength = 4+3*ncolors;
      break;
    }
  case PROP_COMPRESSION: minlength = 1; break;
  case PROP_OPACITY:     minlength = 4; break;
  case PROP_APPLY_MASK:  minlength = 4; break;
  case PROP_OFFSETS:     minlength = 8; break;
  case PROP_MODE:        minlength = 4; break;
  default:               minlength = 0; break;
  }
  if( length < minlength ) {
    FatalBadXCF("Short %s property at %" PRIX32 " (%" PRIu32 "<%" PRIu32 ")",
                showPropType(type),ptr,length,minlength);
    return XCF_ERROR;
  }
  *master = ptr+8+length ;
  total = 8 + length + (type != PROP_END ? 8 : 0) ;
  if( total < length ) { /* Check overwrap */
    FatalBadXCF("Overlong property at %" PRIX32, ptr);
    return XCF_ERROR;
  }
  if((response = xcfCheckspace(ptr,total,"Overlong property at %" PRIX32,ptr)) != 0) {
      return XCF_ERROR;
  }
  *typeOut = type;
  return XCF_OK;
}

const char*
xcfString(uint32_t ptr,uint32_t *after)
{
  uint32_t length ;
  unsigned i ;
  ICONV_CONST char *utf8master ;

  if (xcfCheckspace(ptr,4,"(string length)") != XCF_OK) {
      return XCF_PTR_EMPTY;
  }
  length = xcfL(ptr) ;
  ptr += 4 ;
  if (xcfCheckspace(ptr,length,"(string)") != XCF_OK) {
      return XCF_PTR_EMPTY;
  }
  utf8master = (ICONV_CONST char*)(xcf_file+ptr) ;
  if( after ) *after = ptr + length ;
  if( length == 0 || utf8master[length-1] != 0 ) {
    FatalBadXCF("String at %" PRIX32 " not zero-terminated",ptr-4);
    return XCF_PTR_EMPTY;
  }
  length-- ;

  if( use_utf8 ) return utf8master ;

  /* We assume that the local character set includes ASCII...
   * Check if conversion is needed at all
   */
  for( i=0 ; ; i++ ) {
    if( i == length )
      return utf8master ; /* Only ASCII after all */
    if( utf8master[i] == 0 ) {
      FatalBadXCF("String at %" PRIX32 " has embedded zeroes",ptr-4);
      return XCF_PTR_EMPTY;
    }
    if( (int8_t) utf8master[i] < 0 )
      break ;
  }
#ifdef HAVE_ICONV
  {
    size_t targetsize = length+1 ;
    int sloppy_translation = 0 ;
    iconv_t cd = iconv_open("//TRANSLIT","UTF-8");
    if( cd == (iconv_t) -1 ) {
      cd = iconv_open("","UTF-8");
      sloppy_translation = 1 ;
    }
    if( cd == (iconv_t) -1 )
      iconv_close(cd) ; /* Give up; perhaps iconv doesn't know UTF-8 */
    else
      while(1) {
        char *buffer = xcfmalloc(targetsize) ;
        ICONV_CONST char *inbuf = utf8master ;
        char *outbuf = buffer ;
        size_t incount = length ;
        size_t outcount = targetsize ;
        while(1) { /* Loop for systems without //ICONV support */
          size_t result = iconv(cd,&inbuf,&incount,&outbuf,&outcount) ;
          if( result == (size_t)-1 && errno == EILSEQ &&
              sloppy_translation && outcount > 0 ) {
            *outbuf++ = '?' ;
            outcount-- ;
            while( (int8_t)*inbuf < 0 ) inbuf++, incount-- ;
            continue ;
          }
          if( result != (size_t)-1 ) {
            if( outcount == 0 )
              errno = E2BIG ;
            else {
              *outbuf = 0 ;
              iconv_close(cd) ;
              return buffer ;
            }
          }
          break ;
        }
        if( errno == EILSEQ || errno == EINVAL ) {
          FatalBadXCF("Bad UTF-8 encoding '%s' at %" PRIXPTR,
                      inbuf,(uintptr_t)((inbuf-utf8master)+ptr));
          return XCF_PTR_EMPTY;
        }
        if( errno == E2BIG ) {
          targetsize += 1+incount ;
          xcffree(buffer) ;
          continue ;
        }
        FatalUnexpected("!iconv on layer name at %" PRIX32, ptr);
        return XCF_PTR_EMPTY:
      }
  }
#endif
  {
    static int warned = 0 ;
    if( !warned ) {
      fprintf(stderr,_("Warning: one or more layer names could not be\n"
                       "         translated to the local character set.\n"));
      warned = 1 ;
    }
  }
  return utf8master ;
}

/* ****************************************************************** */

int computeDimensions(struct tileDimensions *d)
{
  // [ CVE-2019-5086 and CVE-2019-5087 ]
  // This part of code is the check to prevent integer overflow, see CVE-2019-5086 and CVE-2019-5087

  if ((d->c.l + d->width) * 4 > INT_MAX) {
    FatalBadXCF(("Width is too large (%d)! Stopping execution...\n"), (d->c.l + d->width));
    return XCF_ERROR;
  }

  if ((d->c.t + d->height) * 4 > INT_MAX) {
    FatalBadXCF(("Height is too large (%d)! Stopping execution...\n"), (d->c.t + d->height));
    return XCF_ERROR;
  }
  // [ CVE-2019-5086 and CVE-2019-5087 ]

  d->c.r = d->c.l + d->width ;
  d->c.b = d->c.t + d->height ;
  d->tilesx = (d->width+TILE_WIDTH-1)/TILE_WIDTH ;
  d->tilesy = (d->height+TILE_HEIGHT-1)/TILE_HEIGHT ;
  d->ntiles = d->tilesx * d->tilesy ;

  return XCF_OK;
}

struct xcfImage XCF ;

int
getBasicXcfInfo(void)
{
  uint32_t ptr, data, layerfile ;
  PropType type ;
  int i, j ;

  int errorStatus;
  uint32_t ptrout;

  if (xcfCheckspace(0,14+7*4,"(very short)") != XCF_OK) {
     return XCF_ERROR;
  }

  if( strcmp((char*)xcf_file,"gimp xcf file") == 0 )
    XCF.version = 0 ;
  else if( xcf_file[13] == 0 &&
          sscanf((char*)xcf_file,"gimp xcf v%d",&XCF.version) == 1 )
    ;
  else {
    FatalBadXCF(_("Not an XCF file at all (magic not recognized)"));
    return XCF_ERROR;
  }

  if (XCF.version < 0 || XCF.version > 2) {
    fprintf(stderr, _("Warning: XCF version %d not supported (trying anyway...)\n"), XCF.version);
  }

  XCF.compression = COMPRESS_NONE ;
  XCF.colormapptr = 0 ;

  ptr = 14 ;
  XCF.width    = xcfL(ptr); ptr += 4 ;
  XCF.height   = xcfL(ptr); ptr += 4 ;
  XCF.type     = xcfL(ptr); ptr += 4 ;
  while( (errorStatus = xcfNextprop(&ptr,&data, &type)) != XCF_ERROR && type != PROP_END ) {
      if (errorStatus != XCF_OK) {
          return XCF_ERROR;
      }

    switch(type) {
    case PROP_COLORMAP:
      XCF.colormapptr = data ;
      break ;
    case PROP_COMPRESSION:
      XCF.compression = xcf_file[data] ;
      break ;
    default:
      /* Ignore unknown properties */
      break ;
    }
  }

  layerfile = ptr ;
  XCF.numLayers = 0;
  while (1) {
      errorStatus = xcfOffset(ptr,8*4, &ptrout);
      if (errorStatus != XCF_OK) {
          return XCF_ERROR;
      }
      if (!ptrout) {
          break;
      }
      XCF.numLayers++;
      ptr+=4;
  }
  XCF.layers = xcfmalloc(XCF.numLayers * sizeof(struct xcfLayer)) ;
  for( i = 0 ; i < XCF.numLayers ; i++ ) {
    struct xcfLayer *L = XCF.layers + i ;
    ptr = xcfL(layerfile+4*(XCF.numLayers-1-i)) ;

    L->mode = GIMP_NORMAL_MODE ;
    L->opacity = 255 ;
    L->isVisible = 1 ;
    L->hasMask = 0 ;
    L->dim.width = xcfL(ptr); ptr+=4 ;
    L->dim.height = xcfL(ptr); ptr+=4 ;
    L->type = xcfL(ptr); ptr+=4 ;
    L->name = xcfString(ptr,&ptr);
    if (L->name == XCF_PTR_EMPTY) {
        return XCF_ERROR;
    }
    L->propptr = ptr ;

    L->isGroup = 0;
    L->pathLength = 0;
    L->path = NULL;

    while( (errorStatus = xcfNextprop(&ptr,&data, &type)) != XCF_ERROR && type != PROP_END ) {
        if (errorStatus != XCF_OK) {
            xcffree(XCF.layers);
            XCF.layers = XCF_PTR_EMPTY;
            return XCF_ERROR;
        }
      switch(type) {
      case PROP_OPACITY:
        L->opacity = xcfL(data);
        if( L->opacity > 255 )
          L->opacity = 255 ;
        break ;
      case PROP_VISIBLE:
        L->isVisible = xcfL(data) != 0 ;
        break ;
      case PROP_APPLY_MASK:
        L->hasMask = xcfL(data) != 0 ;
        break ;
      case PROP_OFFSETS:
        L->dim.c.l = (int32_t)(xcfL(data  )) ;
        L->dim.c.t = (int32_t)(xcfL(data+4)) ;
        break ;
      case PROP_MODE:
        L->mode = xcfL(data);
        break ;
      case PROP_GROUP_ITEM:
        L->isGroup = 1 ;
        break;
      case PROP_ITEM_PATH:
        L->pathLength = (ptr - data - 2) / 4 ;

        if ( L->pathLength != 0 ) {

          L->path = xcfmalloc( L->pathLength * sizeof(unsigned) ) ;

          for ( j = 0; j!=L->pathLength; j++ )
            *(L->path + j) = (unsigned)xcfL(data + 4 * j);
        }
        break;
      default:
        /* Ignore unknown properties */
        break ;
      }
    }
    if ((errorStatus = xcfCheckspace(ptr,8,"(end of layer %s)",L->name)) != XCF_OK) {
        xcffree(XCF.layers);
        XCF.layers = XCF_PTR_EMPTY;
        return XCF_ERROR;
    }
    L->pixels.tileptrs = 0 ;
    if (xcfOffset(ptr  , 4*4, &(L->pixels.hierarchy)) != XCF_OK) {
        xcffree(XCF.layers);
        XCF.layers = XCF_PTR_EMPTY;
        return XCF_ERROR;
    }
    L->mask.tileptrs = 0 ;
    if (xcfOffset(ptr+4, 4*4, &(L->mask.hierarchy)) != XCF_OK) {
        xcffree(XCF.layers);
        XCF.layers = XCF_PTR_EMPTY;
        return XCF_ERROR;
    }

    if (computeDimensions(&L->dim) != XCF_OK) {
        return XCF_ERROR;
    }
  }
  return XCF_OK;
}
