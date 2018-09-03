import ImpulseResponse

using Test
@everywhere using DeviceUnderTest
@everywhere using Soundcard


let (f,h,d,t) = ImpulseResponse.expsinesweep_simulate()
end


function asio_test()
    ms = zeros(1,8)
    ms[1,1] = 1
    mm = zeros(8,1)
    mm[1,1] = 1
    ImpulseResponse.expsinesweep_asio(ms, mm, 48000, 22, 22000, 10, 3, [([1.0],[1.0])], -6)    
end
let (f,h,d,t) = asio_test()
end


function fileio_test()
    ms = zeros(1,2)
    ms[1,1] = 1
    mm = zeros(8,8)
    for i = 1:8
        mm[i,i] = 1
    end
    r = DeviceUnderTest.register()
    expsinesweep_fileio(r[0x9fefe994b7e95bf1], ms, mm, 47999.6, 22, 22000, 10, 3, [([1.0],[1.0])], -6, -10, 800, 2000, 0.5, 2, 3)
end



function fileio_asio_test()
    ms = zeros(1,2)
    ms[1,2] = 1
    mm = zeros(8,1)
    mm[1,1] = 1
    r = DeviceUnderTest.register()
    ImpulseResponse.expsinesweep_fileio_asio(r[0x9fefe994b7e95bf1], ms, mm, 48000, 47999.6, 22, 22000, 10, 3, [([1.0],[1.0])], -6, -10, 800, 2000, 0.5, 2, 3)
end



function asio_fileio_test()
    ms = zeros(1,8)
    ms[1,1] = 1
    mm = zeros(8,8)
    for i = 1:8
        mm[i,i] = 1
    end
    r = DeviceUnderTest.register()
    ImpulseResponse.expsinesweep_asio_fileio(r[0x9fefe994b7e95bf1], ms, mm, 48000, 47999.6, 22, 22000, 10, 3, [([1.0],[1.0])], -6, -10, 800, 2000, 0.5, 2, 3)
end


function clockdrift_test()
    ms = zeros(1,2)
    ms[1,1] = 1
    mm = zeros(8,1)
    mm[1,1] = 1
    r = DeviceUnderTest.register()
    rept = 1
    ImpulseResponse.measureclockdrift(r[0x9fefe994b7e95bf1], ms, mm, rept, 48000, 800, 2000, 0.5, -6)
end