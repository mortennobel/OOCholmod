solution "oocholmod_solution"
  configurations { "Release", "Debug" }
  libdirs { "/usr/local/lib/" }
  location "project" 
  platforms { "Universal" }
  
project "oocholmod"
  kind     "StaticLib"
  language "C++"
  files    { "src/**.h","src/**.cpp"}		
  includedirs { "src", "/usr/local/include"}
  buildoptions "-std=c++11 -stdlib=libc++"
  links { "Accelerate.framework" }
  buildoptions { "-F/Library/Frameworks/" , "-F/usr/local/lib/", "-F./lib/debug"}
  linkoptions { "-F/Library/Frameworks/" , "-F/usr/local/lib/", "-lcholmod", "-lccolamd", "-lcolamd", "-lsuitesparseconfig", "-lcamd","-lamd",  "-lccolamd", "-lamd"}
  
  configuration "Debug"
    targetdir "lib/debug"
    defines { "DEBUG" }
    flags { "Symbols" }		

  configuration "Release"
    targetdir "lib/release"
    defines { "NDEBUG" }
    flags { "Optimize" } 	
	
  if _ACTION == "clean" then
    os.rmdir("lib")
  end
  
project "oocholmod_unittest"
  kind "ConsoleApp"
  language "C++"
  files    {"test/**.h","test/**.cpp"}
  includedirs { "src", "test", "/usr/local/include" }
  buildoptions "-std=c++11 -stdlib=libc++"
  links { "Accelerate.framework"}
  buildoptions { "-F/Library/Frameworks/" , "-F/usr/local/lib/", "-F./lib/debug"}
  linkoptions { "-F/Library/Frameworks/" , "-F/usr/local/lib/", "-lcholmod", "-lccolamd", "-lcolamd", "-lsuitesparseconfig", "-lcamd","-lamd",  "-lccolamd", "-lamd"}
  configuration "Debug"
    targetdir "bin/debug"
    defines { "DEBUG"}
    flags { "Symbols" }	

  configuration "Release"
    targetdir "bin/release"
    defines { "NDEBUG"}
    flags { "Optimize" }
	
  if _ACTION == "clean" then
    os.rmdir("bin")
  end
