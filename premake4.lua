solution "res"
    includedirs { "src","src/**" }
    
    configurations{ "Debug", "Release" }
        configuration "Debug"
            defines { "DEBUG" }
            flags { "Symbols" }

        configuration "Release"
            defines { "NDEBUG" }
            flags { "Optimize" }

    configuration { "linux", "gmake" }
        buildoptions { "-O3", "-pthread", "-Wall", "-Wextra", 
                       "-ansi", "-pedantic", "-std=c++11" }
        links { "pthread" }
        linkoptions { }

    project "res"
        kind "ConsoleApp"
        language "C++"
        files { "src/**.h", "src/**.cpp" }

