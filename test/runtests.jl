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


function minimalphase_eq_test()
    
    ms = [0 0 0 0 1 0 0 0]
    mm = [1,0]
    mm = mm[:,:]
    fs = 48000
    atten = -20
    n = 1024

    # transfer function without any equalization
    f0 = 22
    f1 = 22000
    f, m, ϕ, xds = ImpulseResponse.minimalphase_fir_original(ms, mm, fs, f0, f1, atten, n);
    

    # tuning process
    f0 = 142
    f1 = 16000
    fx = 214
    h = ImpulseResponse.minimalphase_fir_design(f, m, ϕ, xds, fs, f0, f1, fx, n);
    # ...
    f1 = 13564
    h = ImpulseResponse.minimalphase_fir_design(f, m, ϕ, xds, fs, f0, f1, fx, n);


    # plot the verification
    f0 = 22
    f1 = 22000
    ImpulseResponse.minimalphase_fir_verification(h, ms, mm, fs, f0, f1, atten, n);
end