import ImpulseResponse
using Test
using DeviceUnderTest


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


function fileio_asio_test()
    ms = zeros(1,8)
    ms[1,1] = 1
    mm = zeros(8,1)
    mm[1,1] = 1
    r = DeviceUnderTest.register()
    ImpulseResponse.expsinesweep_fileio_asio(r[0x9fefe994b7e95bf1], ms, mm, 48000, 47999.6, 22, 22000, 10, 3, [([1.0],[1.0])], -6, -10, 800, 2000, 0.5, 2, 3)
end