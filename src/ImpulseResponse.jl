module ImpulseResponse
using Statistics
using Random
using LinearAlgebra
using Distributed
using SharedArrays
using Libaudio
using Soundcard
using DeviceUnderTest
using WAV



"""
simulation
"""
function expsinesweep_simulate(fs=48000., f0=22, f1=22000, ts=10, td=3, eq=[([1.0],[1.0])])
    s = Libaudio.expsinesweep(f0, f1, ts, fs)
    m = length(s)
    n = round(Int, td * fs)
    x = zeros(m+n, 1)
    x[1:m,1] = s
    for i in eq
        x = Libaudio.filt(i[1], i[2], x)  
    end
    ellip_b = [0.165069005881145, 0.064728220211450, 0.237771476924023, 0.237771476924022, 0.064728220211450, 0.165069005881146]
    ellip_a = [1.0, -1.544070855644442, 2.377072040756431, -1.638501402700271, 0.950992608325718, -0.210354984704200]
    r = Libaudio.filt(ellip_b, ellip_a, x)
    p = median(abs.(r))
    r[r.>p] .= p
    return Libaudio.impulse(s, n, f0, f1, fs, r)
end



"""
soundcard -> soundcard measurement

# Arguments
- 'ms::Matrix': mixing matrix of loudspeakers
- 'mm::Matrix': mixing matrix of microphones
- 'fs': sample rate of the soundcard
- 'f0': starting frequency
- 'f1': stop frequency
- 'ts': ess signal length in seconds
- 'td': ess decaying length in seconds
- 'eq': filter array for equalization
- 'atten': ess signal attenuation
"""
function expsinesweep_asio(ms::Matrix, mm::Matrix, fs=48000, f0=22, f1=22000, ts=10, td=3, eq=[([1.0],[1.0])], atten=-6)
    @assert f0 < f1
    f0 < 1.0 && (f0 = 1.0)
    f1 > fs/2 && (f1 = fs/2)

    s = Libaudio.expsinesweep(f0, f1, ts, fs)
    m = length(s)
    n = round(Int, td * fs)
    x = zeros(m+n, 1)
    x[1:m,1] = s
    for i in eq
        x = Libaudio.filt(i[1], i[2], x) 
    end
    r = Soundcard.playrecord(x * 10^(atten/20), ms, mm, fs)
    return Libaudio.impulse(s, n, f0, f1, fs, r)
end



"""
dut -> dut measurement

# Arguments
- 'f': dut registered function dictionary
- 'ms::Matrix': mixing matrix of loudspeakers
- 'mm::Matrix': mixing matrix of microphones
- 'fs': sample rate of the dut
- 'f0': starting frequency of the ess signal
- 'f1': stop frequency of the ess signal
- 'ts': ess signal length in seconds
- 'td': ess decaying length in seconds
- 'eq': filter array for equalization
- 'atten': attenuation of the ess signal in dB
- 'syncatten': attenuation of the sync symbol in dB
- 'f2': starting frequency of the sync symbol
- 'f3': stop frequency of the sync symbol
- 'tsm': time for syncsymbol generation
- 'tsd': time for syncsymbol decaying
- 'tcs': time for context switching

# Examples
    using DeviceUnderTest
    f = DeviceUnderTest.register()
    expsinesweep_fileio(f[some device code], ...)
"""
function expsinesweep_fileio(f, ms::Matrix, mm::Matrix, fs=47999.6, f0=22, f1=22000, ts=10, td=3, eq=[([1.0],[1.0])], atten=-6, syncatten=-10, f2=800, f3=2000, tsm=0.5, tsd=2, tcs=3)
    @assert f0 < f1
    f0 < 1.0 && (f0 = 1.0)
    f1 > fs/2 && (f1 = fs/2)

    s = Libaudio.expsinesweep(f0, f1, ts, fs)
    m = length(s)
    n = round(Int, td * fs)
    x = zeros(m+n, 1)
    x[1:m,1] = s
    for i in eq
        x = Libaudio.filt(i[1], i[2], x)  
    end
    sym = convert(Vector{Float64}, Libaudio.symbol_expsinesweep(f2, f3, tsm, fs))
    y = Libaudio.encode_syncsymbol(tcs, sym, tsd, x * 10^(atten/20), fs, 1, syncatten)
    out = randstring() * ".wav"
    av = [8000, 16000, 44100, 48000, 96000, 192000]
    wavwrite(DeviceUnderTest.mixer(y, ms), out, Fs=av[findmin(abs.(av.-fs))[2]], nbits=32)

    try
        f[:init]()
        f[:readyplayrecord](out, ceil(size(y,1)/fs), true)
        f[:playrecord](true)
        DeviceUnderTest.raw2wav16bit("./capture/mic_pcm_before_resample_8ch_48000.raw", size(mm,1), 48000, "./capture/mic_pcm_before_resample_8ch_48000.wav")
        r, rate = wavread("./capture/mic_pcm_before_resample_8ch_48000.wav")

        nx = size(x,1)
        p = Libaudio.decode_syncsymbol(r, sym, tsd, nx/fs, fs)
        c = size(r,2)
        rd = zeros(nx, c)
        tol = 2048
        for i = 1:c
            l = p[i] - tol
            rd[:,i] = r[l:l+nx-1,i]
        end
        return Libaudio.impulse(s, n, f0, f1, fs, rd)
    finally
        rm(out, force=true)
        rm("./capture/", force=true, recursive=true)
    end
end


"""
soundcard -> dut measurement
"""
function expsinesweep_asio_fileio(f, ms::Matrix, mm::Matrix, fs=48000, fm=47999.6, f0=22, f1=22000, ts=10, td=3, eq=[([1.0],[1.0])], atten=-6, syncatten=-10, f2=800, f3=2000, tsm=0.5, tsd=2, tcs=3)
    @assert nprocs() > 1
    @assert f0 < f1
    f0 < 1.0 && (f0 = 1.0)
    f1 > fs/2 && (f1 = fs/2)
    wid = workers()

    s = Libaudio.expsinesweep(f0, f1, ts, fs)
    m = length(s)
    n = round(Int, td * fs)
    x = zeros(m+n, 1)
    x[1:m,1] = s
    for i in eq
        x = Libaudio.filt(i[1], i[2], x)  
    end
    sym = convert(Vector{Float64}, Libaudio.symbol_expsinesweep(f2, f3, tsm, fs))
    y = Libaudio.encode_syncsymbol(tcs, sym, tsd, x * 10^(atten/20), fs, 1, syncatten)

    try
        f[:init]()
        f[:readyrecord](ceil(size(y,1)/fs), true)
        phy = Soundcard.mixer(y, ms)
        pcm = SharedArray{Float32,1}(Soundcard.to_interleave(phy))
        # done = remotecall(Soundcard.play, wid[1], y, ms, fs) # latency is high
        done = remotecall(Soundcard.play, wid[1], size(phy), pcm, fs)  # latency is low
        f[:record](true)
        fetch(done)
        DeviceUnderTest.raw2wav16bit("./capture/mic_pcm_before_resample_8ch_48000.raw", size(mm,1), 48000, "./capture/mic_pcm_before_resample_8ch_48000.wav")
        r,rate = wavread("./capture/mic_pcm_before_resample_8ch_48000.wav")
        
        # decode async signal
        nx = size(x,1)
        tx = nx / fs
        syma = convert(Vector{Float64}, Libaudio.symbol_expsinesweep(f2, f3, tsm, fm))
        p = Libaudio.decode_syncsymbol(r, syma, tsd, tx, fm)
        sa = Libaudio.expsinesweep(f0, f1, ts, fm)
        nxa = round(Int, tx * fm)
        na = nxa-length(sa)
        @info "decay samples derived from (td x fm) and (tx x fm)" round(Int,td * fm) na
        c = size(r,2)
        tol = 2048
        rd = zeros(nxa, c)
        for i = 1:c
            l = p[i] - tol
            rd[:,i] = r[l:l+nxa-1,i]
        end
        return Libaudio.impulse(sa, na, f0, f1, fm, rd)
    finally
        rm("./capture/", force=true, recursive=true)
    end  
end


"""
dut -> soundcard measurement
"""
function expsinesweep_fileio_asio(f, ms::Matrix, mm::Matrix, fs=48000, fm=47999.6, f0=22, f1=22000, ts=10, td=3, eq=[([1.0],[1.0])], atten=-6, syncatten=-10, f2=800, f3=2000, tsm=0.5, tsd=2, tcs=3)
    @assert nprocs() > 1
    @assert f0 < f1
    f0 < 1.0 && (f0 = 1.0)
    f1 > fs/2 && (f1 = fs/2)
    wid = workers()

    s = Libaudio.expsinesweep(f0, f1, ts, fm)
    m = length(s)
    n = round(Int, td * fm)
    x = zeros(m+n, 1)
    x[1:m,1] = s
    for i in eq
        x = Libaudio.filt(i[1], i[2], x)  
    end
    sym = convert(Vector{Float64}, Libaudio.symbol_expsinesweep(f2, f3, tsm, fm))
    y = Libaudio.encode_syncsymbol(tcs, sym, tsd, x * 10^(atten/20), fm, 1, syncatten)
    out = randstring() * ".wav"
    wavwrite(DeviceUnderTest.mixer(y, ms), out, Fs=fs, nbits=32)
    
    try
        f[:init]()
        f[:readyplay](out)
        done = remotecall(f[:play], wid[1])
        r = Soundcard.record(round(Int,fs*size(y,1)/fm), mm, fs)
        fetch(done)

        nx = size(x,1)
        tx = nx / fm
        syma = convert(Vector{Float64}, Libaudio.symbol_expsinesweep(f2, f3, tsm, fs))
        p = Libaudio.decode_syncsymbol(r, syma, tsd, tx, fs)
        sa = Libaudio.expsinesweep(f0, f1, ts, fs)
        nxa = round(Int, tx * fs)
        na = nxa-length(sa)
        @info "decay samples derived from (td x fs) and (tx x fs)" round(Int,td * fs) na
        c = size(r,2)
        tol = 2048
        rd = zeros(nxa, c)
        for i = 1:c
            l = p[i] - tol
            rd[:,i] = r[l:l+nxa-1,i]
        end
        return Libaudio.impulse(sa, na, f0, f1, fs, rd)
    finally
        rm(out, force=true)
    end
end




# const sndcard_n_out = 8
# const sndcard_n_in = 8
# const mic_n = 1
# mixplay = zeros(Float32, size(essd,2), sndcard_n_out)
# mixplay[1,2] = 1.0f0
# mixrec = zeros(Float32, sndcard_n_in, mic_n)
# mixrec[2,1] = 1.0f0
#
# note: adjust atten for the level of the ess signal
# note: adjust syncatten for the level of the sync symbol, if too high mixer would prompt with mic clipping error!
# note: use long t_ess when possible, for the build-up of the stimulus energy
# note: use long t_decay if room or system dunamics are reverberant
# note: eq = [(b,a),(b,a)...] is prepending transfer function for filter verification, must be designed according to fs or fsd based on mode
# note: mode[1] is the physical device for playback, mode[2] is the physical device for recording






# mixspk = zeros(1,2)
# mixspk[1,1] = 1.0
# mixmic = zeros(9,1)
# mixmic[9,1] = 1.0
function measureclockdrift(f, ms::Matrix{Float64}, mm::Matrix{Float64}, rep=3, fs=48000, f2=800, f3=2000, tsm=0.5, atten=-6)
    # 
    #   +--------+--------+--------+--------+--------+ => 5 samples in digital domain played via dut's speaker, whose sample interval is Td.
    #   +-----+-----+-----+-----+-----+-----+-----+--- => 7 samples captured by the standard sampler of the soundcard,
    #  
    #   5 x Td ≈ 7 x Tr
    #   or formly, N/Fd ≈ Nm / Fr
    #   Fd ≈ N / Nm x Fr
    #   Fd/Fr ≈ N/Nm = 5/7  
    #
    @assert nprocs() > 1
    wpid = workers()
    @info "start measure clock drift"

    # fileio -> asio
    sync = 10^(atten/20) * convert(Vector{Float64}, Libaudio.symbol_expsinesweep(f2, f3, tsm, fs))
    @info "sync samples" length(sync)
    period = [zeros(round(Int,100fs),1); sync]
    signal = [zeros(round(Int,3fs),1); sync; repeat(period,rep,1); zeros(round(Int,3fs),1)]
    @info "signal train formed"

    f[:init]()
    out = randstring() * ".wav"
    playback = "dutplaybackclockdriftmeasure.wav"
    wavwrite(Device.mixer(signal, ms), out, Fs=fs, nbits=32)
    @info "filesize in MiB" filesize(out)/1024/1024
    f[:readyplay](out)
    @info "singal pushed to device"

    done = remotecall(f[:play], wpid[1])
    r = Soundcard.record(size(signal,1), mm, fs)
    fetch(done)
    wavwrite(r, "clockdrift.wav", Fs=fs, nbits=32)
    @info "recording written to clockdrift.wav"

    # syncs = 10^(-6/20) * LibAudio.syncsymbol(800, 2000, 1, fss)
    # info("  syncs samples: $(length(syncs))")
    # note: sync is approximately invariant due to its short length
    measure = Array{Tuple{Float64, Float64, Float64},1}()
    for k = 1:size(r,2)
        lbs,pk,pkf,y = Libaudio.extractsymbol(r[:,k], sync, rep+1)
        pkfd = diff(pkf)
        chrodrift_100sec = ((pkfd[end] - pkfd[1]) / (rep-1))/fs
        freqdrift_100sec = (size(period,1) - median(pkfd))/fs
        push!(measure, (fs * size(period,1)/median(pkfd), freqdrift_100sec, chrodrift_100sec))
    end
    measure
end

end # module
