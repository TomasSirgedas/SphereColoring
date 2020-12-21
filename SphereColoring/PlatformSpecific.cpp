#include "PlatformSpecific.h"
#include <Windows.h>


bool isKeyDown( char c )
{
   return GetKeyState( (int) c ) < 0;
}
