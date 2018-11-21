module ImpulseResponse

using Statistics
using Random
using LinearAlgebra
using Distributed
using SharedArrays

using PyPlot

using Libaudio
using Soundcard
using DeviceUnderTest




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
    Libaudio.wavwrite(out, DeviceUnderTest.mixer(y, ms), av[findmin(abs.(av.-fs))[2]], 32)

    root = joinpath(Libaudio.folder(), Libaudio.logfile())
    try
        f[:init]()
        f[:readyplayrecord](out, ceil(size(y,1)/fs), true)
        f[:playrecord](true)
        r, rate = Libaudio.wavread("./raw_out_mic_all_16bit_48k_8ch_mic_pcm_before_resample_8ch_48000.wav", Float64)
        nx = size(x,1)
        val, p = Libaudio.decode_syncsymbol(r, sym, tsd, nx/fs, fs)
        val || Libaudio.printl(root, :light_red, Libaudio.nows() * " | ImpulseResponse.expsinesweep_fileio: decode symbol failure")
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
        rm("./raw_out_mic_all_16bit_48k_8ch_mic_pcm_before_resample_8ch_48000.wav", force=true)
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

    root = joinpath(Libaudio.folder(), Libaudio.logfile())
    try
        f[:init]()
        f[:readyrecord](ceil(size(y,1)/fs), true)
        phy = Soundcard.mixer(y, ms)
        pcm = SharedArray{Float32,1}(Soundcard.interleave(phy))
        done = remotecall(Soundcard.play, wid[1], size(phy), pcm, fs)  # latency is low
        f[:record](true)
        fetch(done)
        r,rate = Libaudio.wavread("./raw_out_mic_all_16bit_48k_8ch_mic_pcm_before_resample_8ch_48000.wav", Float64)
        nx = size(x,1)
        tx = nx / fs
        syma = convert(Vector{Float64}, Libaudio.symbol_expsinesweep(f2, f3, tsm, fm))
        val, p = Libaudio.decode_syncsymbol(r, syma, tsd, tx, fm)
        val || Libaudio.printl(root, :light_red, Libaudio.nows() * " | ImpulseResponse.expsinesweep_asio_fileio: decode symbol failure")
        sa = Libaudio.expsinesweep(f0, f1, ts, fm)
        nxa = round(Int, tx * fm)
        na = nxa-length(sa)

        Libaudio.printl(root, :light_cyan, Libaudio.nows() * " | ImpulseResponse.expsinesweep_asio_fileio: decay samples [td x fm]/[tx x fm] $(round(Int,td * fm))/$(na)")
        c = size(r,2)
        tol = 2048
        rd = zeros(nxa, c)
        for i = 1:c
            l = p[i] - tol
            rd[:,i] = r[l:l+nxa-1,i]
        end
        return Libaudio.impulse(sa, na, f0, f1, fm, rd)
    finally
        rm("./raw_out_mic_all_16bit_48k_8ch_mic_pcm_before_resample_8ch_48000.wav", force=true)
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
    Libaudio.wavwrite(out, DeviceUnderTest.mixer(y, ms), fs, 32)
    
    root = joinpath(Libaudio.folder(), Libaudio.logfile())
    try
        f[:init]()
        f[:readyplay](out)
        done = remotecall(f[:play], wid[1])
        r = convert(Matrix{Float64}, Soundcard.record(round(Int,fs*size(y,1)/fm), mm, fs))
        fetch(done)
        nx = size(x,1)
        tx = nx / fm
        syma = convert(Vector{Float64}, Libaudio.symbol_expsinesweep(f2, f3, tsm, fs))
        val, p = Libaudio.decode_syncsymbol(r, syma, tsd, tx, fs)
        val || Libaudio.printl(root, :light_red, Libaudio.nows() * " | ImpulseResponse.expsinesweep_fileio_asio: decode symbol failure")
        sa = Libaudio.expsinesweep(f0, f1, ts, fs)
        nxa = round(Int, tx * fs)
        na = nxa-length(sa)

        Libaudio.printl(root, :light_cyan, Libaudio.nows() * " | ImpulseResponse.expsinesweep_fileio_asio: decay samples [td x fs]/[tx x fs] $(round(Int,td * fs))/$(na)")        
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
"""
dut@fm -> soundcard@fs
"""
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

    root = joinpath(Libaudio.folder(), Libaudio.logfile())
    Libaudio.printl(root, :light_cyan, Libaudio.nows() * " | ImpulseResponse.measureclockdrift: start measuring device clock drift")

    # fileio -> asio
    sync = 10^(atten/20) * convert(Vector{Float64}, Libaudio.symbol_expsinesweep(f2, f3, tsm, fs))
    Libaudio.printl(root, :light_cyan, Libaudio.nows() * " | ImpulseResponse.measureclockdrift: measureclockdrift: sync samples $(length(sync))")

    period = [zeros(round(Int,100fs),1); sync]
    signal = [zeros(round(Int,3fs),1); sync; repeat(period,rep,1); zeros(round(Int,3fs),1)]
    Libaudio.printl(root, :light_cyan, Libaudio.nows() * " | ImpulseResponse.measureclockdrift: stimulus formed")

    out = randstring() * ".wav"
    Libaudio.wavwrite(out, DeviceUnderTest.mixer(signal, ms), fs, 32)
    Libaudio.printl(root, :light_cyan, Libaudio.nows() * " | ImpulseResponse.measureclockdrift: filesize $(filesize(out)/1024/1024) MiB")

    measure = Array{Tuple{Float64, Float64, Float64},1}()
    try
        f[:init]()
        f[:readyplay](out)
        Libaudio.printl(root, :light_cyan, Libaudio.nows() * " | ImpulseResponse.measureclockdrift: stimulus pushed to device")

        done = remotecall(f[:play], wpid[1])
        r = convert(Matrix{Float64}, Soundcard.record(size(signal,1), mm, fs))
        fetch(done)
        # Libaudio.wavwrite("clockdrift.wav", r, fs, 16)
        val = true
        for k = 1:size(r,2)
            flag,lbs,pk,pkf,y = Libaudio.extractsymbol(r[:,k], sync, rep+1)
            val = val && flag
            pkfd = diff(pkf)
            chrodrift_100sec = ((pkfd[end] - pkfd[1]) / (rep-1))/fs
            freqdrift_100sec = (size(period,1) - median(pkfd))/fs
            push!(measure, (fs * size(period,1)/median(pkfd), freqdrift_100sec, chrodrift_100sec))
        end
        Libaudio.printl(root, :light_cyan, Libaudio.nows() * " | ImpulseResponse.measureclockdrift: status $val")
    finally
        rm(out)
    end

    for k in measure
        Libaudio.printl(root, :light_cyan, Libaudio.nows() * " | ImpulseResponse.measureclockdrift: time drift every 100 seconds $(k[2]/100)")
        Libaudio.printl(root, :light_cyan, Libaudio.nows() * " | ImpulseResponse.measureclockdrift: temperature drift every 100 seconds $(k[3]/100)")
        Libaudio.printl(root, :light_cyan, Libaudio.nows() * " | ImpulseResponse.measureclockdrift: estimated sample rate $(k[1])")
    end

    measure
end



"""
soundcard@fs -> dut@fm
"""
function measureclockdrift2(f, ms::Matrix{Float64}, mm::Matrix{Float64}, rep=3, fs=48000, f2=800, f3=2000, tsm=0.5, atten=-6)
    # 
    #   +--------+--------+--------+--------+--------+ => 5 samples in digital domain played via standard sampler of the soundcard, whose sample interval is Tr.
    #   +-----+-----+-----+-----+-----+-----+-----+--- => 7 samples captured by the imperfect sampler of the device under test,
    #  
    #   5 x Tr ≈ 7 x Td
    #   or formly, N/Fr ≈ Nm / Fd
    #   Fd ≈ Nm / N x Fr
    #   Fd/Fr ≈ Nm/N = 5/7  
    #
    @assert nprocs() > 1
    wpid = workers()
    root = joinpath(Libaudio.folder(), Libaudio.logfile())


    Libaudio.printl(root, :light_cyan, Libaudio.nows() * " | ImpulseResponse.measureclockdrift2: start measuring device clock drift")

    sync = 10^(atten/20) * convert(Vector{Float64}, Libaudio.symbol_expsinesweep(f2, f3, tsm, fs))
    Libaudio.printl(root, :light_cyan, Libaudio.nows() * " | ImpulseResponse.measureclockdrift2: sync samples $(length(sync))")

    period = [zeros(round(Int,100fs),1); sync]
    signal = [zeros(round(Int,3fs),1); sync; repeat(period,rep,1); zeros(round(Int,3fs),1)]
    Libaudio.printl(root, :light_cyan, Libaudio.nows() * " | ImpulseResponse.measureclockdrift2: stimulus formed")

    measure = Array{Tuple{Float64, Float64, Float64},1}()
    try
        f[:init]()
        f[:readyrecord](ceil(size(signal,1)/fs), true)
        phy = Soundcard.mixer(signal, ms)
        pcm = SharedArray{Float32,1}(Soundcard.interleave(phy))
        done = remotecall(Soundcard.play, wpid[1], size(phy), pcm, fs)  # latency is low
        f[:record](true)
        fetch(done)
        r, rate = Libaudio.wavread("./raw_out_mic_all_16bit_48k_8ch_mic_pcm_before_resample_8ch_48000.wav", Float64)

        val = true
        for k = 1:size(r,2)
            flag,lbs,pk,pkf,y = Libaudio.extractsymbol(convert(Vector{Float64},r[:,k]), sync, rep+1)
            val = val && flag
            pkfd = diff(pkf)
            chrodrift_100sec = ((pkfd[end] - pkfd[1]) / (rep-1))/fs
            freqdrift_100sec = (size(period,1) - median(pkfd))/fs
            push!(measure, (fs * median(pkfd)/size(period,1), freqdrift_100sec, chrodrift_100sec))
        end
        Libaudio.printl(root, :light_cyan, Libaudio.nows() * " | ImpulseResponse.measureclockdrift: status $val")
    finally
        rm("./raw_out_mic_all_16bit_48k_8ch_mic_pcm_before_resample_8ch_48000.wav")
    end
    measure
end






"""
function hFIR = eq_calibration(mix_spk, mix_mic, f_anchor, f_start, f_stop, atten)

    fs = 48000;
    ess_f0 = 10;
    ess_f1 = fs/2;
    ess_time = 30;
    ess_decay = 5;
    %atten = -20;
    nfft = 65536;
    nnyq = nfft/2+1;
    
    %@ measure the un-eq'd impulse response in the anechoic
    [fundamental, harmonics, response_t] = impulse_response_exponential_sine_sweep(mix_spk, mix_mic, ess_f0, ess_f1, ess_time, ess_decay, 'asio', [1], [1], atten);
    
    
    %@ cut the window with the impulse response out for analysis
    p = 4.3e4;
    x = fundamental(p+1:p+nfft);
    x_spec = fft(x)./nfft;
    x_phase = angle(x_spec);
    f = ((0:nnyq-1)'./nfft).*fs;
    x_spec_db = 20*log10(abs(x_spec(1:nnyq)));
    
    % smooth and plots
    Noct = 3;
    x_spec_db_sm = smoothSpectrum(x_spec_db,f,Noct);
    figure(1);
    subplot(2,1,1); plot(20*log10(abs(x))); title('impulse response in time domain'); xlabel('samples'); hold on;
    %subplot(2,1,2); semilogx(f,x_spec_db,f,x_spec_db_sm); title('impulse response in frequency domain, with 1/3 octave smoothing'); xlabel('Hz'); hold on;
    subplot(2,1,2); semilogx(f,x_spec_db_sm); title('impulse response in frequency domain, with 1/3 octave smoothing'); xlabel('Hz'); hold on;
    
    
    % construct a minimal phase filter of the smoothed impulse response
    x_spec_sm = 10.^(x_spec_db_sm/20);
    x_spec_sm = [x_spec_sm; flipud(x_spec_sm(2:end-1))];
    x_spec_sm = x_spec_sm .* exp(x_phase*1i);
    x_sm = real(ifft(x_spec_sm));
    
    zs = fft(x_sm);
    figure(2); semilogx(f, 20*log10(abs(zs(1:nnyq))), 'k'); hold on;
    
    zms = mps(fft(x_sm));
    zm = real(ifft(zms)); % it is not symmetrical
    zms = fft(zm);
    semilogx(f, 20*log10(abs(zms(1:nnyq))), 'c--'); 
    xlabel('Hz'); title('minimal phase filter constructed based on smoothed frequency response'); grid on;
    
    
    % do the actual flatten work
    %f1 = 22;
    f1 = f_start;
    f2 = f_stop;
    f1 = ceil(f1 * nfft / fs);
    f2 = floor(f2 * nfft / fs);
    target = x_spec_db_sm(round(f_anchor * nfft / fs));
    %target = x_spec_db(round(f_anchor * nfft / fs));
    H = zeros(nnyq,1);
    H(f1:f2) = target - x_spec_db_sm(f1:f2);
    %H(f1:f2) = target - x_spec_db(f1:f2);
    H(1:f1-1) = H(f1);
    figure(3); semilogx(f, H, 'r'); hold on; 
    
    % compensation filter in time and frequency domain
    H = 10.^(H/20) * nfft;
    H = [H; flipud(H(2:end-1))];
    H = H .* exp(x_phase*1i);
    h = real(ifft(H));
    hs = fft(h)./nfft;
    semilogx(f, 20*log10(abs(hs(1:nnyq))), 'b--');
    
    hms = mps(fft(h));
    hm = real(ifft(hms)); % it is not symmetrical
    hms = fft(hm)./nfft;
    semilogx(f, 20*log10(abs(hms(1:nnyq))), 'c--'); 
    xlabel('Hz'); title('Compensation filter and its minimal-phase realization'); grid on;
    figure(2); semilogx(f, 20*log10(abs(hms(1:nnyq)))+x_spec_db_sm, 'b'); grid on;
    
    [fundamental, harmonics, response_t] = impulse_response_exponential_sine_sweep(mix_spk, mix_mic, ess_f0, ess_f1, ess_time, ess_decay, 'asio', hm/nfft, [1], atten);
    
    p = 4.3e4;
    x = fundamental(p+1:p+nfft);
    x_spec = fft(x)./nfft;
    x_phase = angle(x_spec);
    f = ((0:nnyq-1)'./nfft).*fs;
    x_spec_db = 20*log10(abs(x_spec(1:nnyq)));
    
    % smooth and plots
    Noct = 3;
    x_spec_db_sm = smoothSpectrum(x_spec_db,f,Noct);
    figure(1);
    subplot(2,1,1); plot(20*log10(abs(x)), 'r'); grid on;
    %subplot(2,1,2); semilogx(f,x_spec_db,f,x_spec_db_sm); grid on;
    subplot(2,1,2); semilogx(f,x_spec_db_sm); grid on;
    
    hFIR = hm / nfft;
end
"""
function minimalphase_eq(ms::Matrix, mm::Matrix, fs, f0, f1=fs/2, fx=150, attenuate=-50, unitcircle=1024, tess=30, tdecay=3, noct=3)
    
    PyPlot.ion()

    fundamental, harmonic, dirac, measure = expsinesweep_asio(ms, mm, fs, f0, f1, tess, tdecay, [([1.0],[1.0])], attenuate)
    f, x = freqz(fundamental[:,1], fs, unitcircle)
    m = length(f)
    ϕ = angle.(x)
    xd = 20log10.(abs.(x))
    xds = Libaudio.smoothspectrum(xd, f, noct)
    
    PyPlot.figure(1)
    PyPlot.semilogx(f, xd)
    PyPlot.semilogx(f, xds)
    PyPlot.title('impulse response in frequency domain, with 1/3 octave smoothing')
    PyPlot.xlabel('Hz')
    PyPlot.grid()

    # construct a minimal phase filter of the smoothed impulse response
    xs = 10.^(xds/20)
    xs = [xs; reverse(xs[2:end-1])]
    fundamental_sm = real(ifft(xs .* exp.(ϕ*im)))
    
    y = fft(fundamental_sm)
    z = Libaudio.mps(y)
    fundamental_mp = real(ifft(z))    # not symmetrical

    figure(2)
    PyPlot.semilogx(f, 20log10.(abs.(y[1:m])))
    PyPlot.semilogx(f, 20log10.(abs.(fft(fundamental_mp)[1:m]))) 
    PyPlot.xlabel('Hz')
    PyPlot.title('minimal phase filter constructed based on smoothed frequency response')
    PyPlot.grid()
    

    # do the actual flatten work
    nfft = 2unitcircle
    ω0 = ceil(f0 * nfft / fs)
    ω1 = floor(f1 * nfft / fs)
    ωx = round(fx * nfft / fs)
    anchor = xds[ωx]

    
    H = zeros(m)
    H[ω0:ω1] .= anchor - xds[ω0:ω1]
    H[1:ω0-1] = H[ω0]
    figure(3) 
    PyPlot.semilogx(f, H)
    
    H = 10.^(H/20) * nfft
    H = [H; reverse(H[2:end-1])]
    H = H .* exp.(ϕ*im)
    h = real(ifft(H))
    hs = fft(h)/nfft
    PyPlot.semilogx(f, 20log10.(abs.(hs[1:m])))
    
    hms = mps(fft(h))
    hm = real(ifft(hms))    # not symmetrical
    hms = fft(hm)/nfft
    PyPlot.semilogx(f, 20log10.(abs.(hms[1:m]))) 
    PyPlot.xlabel('Hz') 
    PyPlot.title('Compensation filter and its minimal-phase realization')
    PyPlot.grid()
    
    figure(2) 
    PyPlot.semilogx(f, 20log10.(abs.(hms[1:m]))+xds) 
    
    hfir = hm / nfft
    fundamental, harmonic, dirac, measure = expsinesweep_asio(ms, mm, fs, f0, f1, tess, tdecay, [(hfir,[1.0])], attenuate)
    f, y = freqz(fundamental[:,1], fs, unitcircle)
    ϕ = angle.(y)
    yd = 20log10.(abs.(y))
    yds = Libaudio.smoothspectrum(yd, f, noct)
    
    PyPlot.figure(1)
    PyPlot.semilogx(f, yd)
    PyPlot.semilogx(f, yds)
    
    return hfir
end


end # module
