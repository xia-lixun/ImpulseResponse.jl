using PyCall
using HDF5
using ImpulseResponse



@pyimport tkinter as tk


speakermix2mat(x) = permutedims([parse(Int,i) for i in split(x[2:end-1])])
microphonemix2mat(x) = [parse(Int,i) for i in split(x[2:end-1],',')]


f = Any
m = Any
ϕ = Any
xds = Any
h = Any      


function measureclicked()
end
function calculateclicked()
end
function validateclicked()
end
function saveclicked()
end


@pydef mutable struct SampleApp <: tk.Tk
    __init__(self, args...; kwargs...) = begin
        tk.Tk[:__init__](self, args...; kwargs...)
        
        
        # label frame parameters
        self[:lf0] = tk.LabelFrame(text="parameters")
        self[:lf0][:pack](fill="both", expand="yes")

        self[:lab00] = tk.Label(self[:lf0], text="sample rate")
        self[:lab00][:pack](side="left")
        self[:ent00v] = tk.StringVar(value="48000")
        self[:ent00] = tk.Entry(self[:lf0], width=10, textvariable=self[:ent00v])
        self[:ent00][:pack](side="left")

        self[:lab01] = tk.Label(self[:lf0], text="nfft/2")
        self[:lab01][:pack](side="left")
        self[:ent01v] = tk.StringVar(value="1024") 
        self[:ent01] = tk.Entry(self[:lf0], width=10, textvariable=self[:ent01v])
        self[:ent01][:pack](side="left")

        self[:lab02] = tk.Label(self[:lf0], text="attenuate(dB)")
        self[:lab02][:pack](side="left")
        self[:ent02v] = tk.StringVar(value="-20")
        self[:ent02] = tk.Entry(self[:lf0], width=10, textvariable=self[:ent02v])
        self[:ent02][:pack](side="left")

        self[:lab03] = tk.Label(self[:lf0], text="speaker mix")
        self[:lab03][:pack](side="left")
        self[:ent03v] = tk.StringVar(value="[0 0 0 0 1 0 0 0]")
        self[:ent03] = tk.Entry(self[:lf0], width=20, textvariable=self[:ent03v])
        self[:ent03][:pack](side="left")

        self[:lab04] = tk.Label(self[:lf0], text="microphone mix")
        self[:lab04][:pack](side="left")
        self[:ent04v] = tk.StringVar(value="[1,0,0,0,0,0,0,0]")
        self[:ent04] = tk.Entry(self[:lf0], width=20, textvariable=self[:ent04v])
        self[:ent04][:pack](side="left")


        # label farme impulse response
        self[:lf1] = tk.LabelFrame(text="impulse response of system WITHOUT equalization")
        self[:lf1][:pack](fill="both", expand="yes")
        self[:lf2] = tk.LabelFrame(text="minimal phase fir equalization filter settings")
        self[:lf2][:pack](fill="both", expand="yes")
        self[:lf3] = tk.LabelFrame(text="impulse response of system WITH equalization")
        self[:lf3][:pack](fill="both", expand="yes")


        self[:lab1] = tk.Label(self[:lf1], text="start(Hz)")
        self[:lab1][:pack](side="left")
        self[:ent1v] = tk.StringVar(value="22")
        self[:ent1] = tk.Entry(self[:lf1], width=10, textvariable=self[:ent1v])
        self[:ent1][:pack](side="left")

        self[:lab2] = tk.Label(self[:lf1], text="stop(Hz)")
        self[:lab2][:pack](side="left")
        self[:ent2v] = tk.StringVar(value="22000")
        self[:ent2] = tk.Entry(self[:lf1], width=10, textvariable=self[:ent2v])
        self[:ent2][:pack](side="left")
 

        self[:btn1] = tk.Button(self[:lf1], text="measure", command=measureclicked)
        self[:btn1][:pack](side="left")


        # label frame filter design
        self[:lab3] = tk.Label(self[:lf2], text="start(Hz)")
        self[:lab3][:pack](side="left")
        self[:ent3v] = tk.StringVar(value="142")
        self[:ent3] = tk.Entry(self[:lf2], width=10, textvariable=self[:ent3v])
        self[:ent3][:pack](side="left")

        self[:lab4] = tk.Label(self[:lf2], text="anchor(Hz)")
        self[:lab4][:pack](side="left")
        self[:ent4v] = tk.StringVar(value="214")
        self[:ent4] = tk.Entry(self[:lf2], width=10, textvariable=self[:ent4v])
        self[:ent4][:pack](side="left")

        self[:lab5] = tk.Label(self[:lf2], text="stop(Hz)")
        self[:lab5][:pack](side="left")
        self[:ent5v] = tk.StringVar(value="13564")
        self[:ent5] = tk.Entry(self[:lf2], width=10, textvariable=self[:ent5v])
        self[:ent5][:pack](side="left")

        self[:btn2] = tk.Button(self[:lf2], text="calculate", command=calculateclicked)
        self[:btn2][:pack](side="left")


        # label frame filter verification
        self[:lab6] = tk.Label(self[:lf3], text="start(Hz)")
        self[:lab6][:pack](side="left")
        self[:ent6v] = tk.StringVar(value="22")
        self[:ent6] = tk.Entry(self[:lf3], width=10, textvariable=self[:ent6v])
        self[:ent6][:pack](side="left")

        self[:lab7] = tk.Label(self[:lf3], text="stop(Hz)")
        self[:lab7][:pack](side="left")
        self[:ent7v] = tk.StringVar(value="22000")
        self[:ent7] = tk.Entry(self[:lf3], width=10, textvariable=self[:ent7v])
        self[:ent7][:pack](side="left")

        self[:btn3] = tk.Button(self[:lf3], text="validate", command=validateclicked)
        self[:btn3][:pack](side="left")


        # label frame save to file
        self[:lf4] = tk.LabelFrame(text="save filter to file")
        self[:lf4][:pack](fill="both", expand="yes")

        self[:lab8] = tk.Label(self[:lf4], text="file name")
        self[:lab8][:pack](side="left")
        self[:ent8v] = tk.StringVar(value="euqalization-20181122.h5")
        self[:ent8] = tk.Entry(self[:lf4], width=40, textvariable=self[:ent8v])
        self[:ent8][:pack](side="left")

        self[:lab9] = tk.Label(self[:lf4], text="dataset name")
        self[:lab9][:pack](side="left")
        self[:ent9v] = tk.StringVar(value="mth1")
        self[:ent9] = tk.Entry(self[:lf4], width=20, textvariable=self[:ent9v])
        self[:ent9][:pack](side="left")

        self[:btn4] = tk.Button(self[:lf4], text="save", command=saveclicked)
        self[:btn4][:pack](side="left")

    end
end
app = SampleApp()


function measureclicked()
    ms = speakermix2mat(app[:ent03][:get]())
    mm = microphonemix2mat(app[:ent04][:get]())
    mm = mm[:,:]
    fs = parse(Int, app[:ent00][:get]())
    n = parse(Int, app[:ent01][:get]())
    atten = parse(Float64, app[:ent02][:get]())

    f0 = parse(Float64, app[:ent1][:get]())
    f1 = parse(Float64, app[:ent2][:get]())

    @info "measurement:" ms mm fs n atten f0 f1
    f_, m_, ϕ_, xds_ = ImpulseResponse.minimalphase_fir_original(ms, mm, fs, f0, f1, atten, n)
    global f = f_
    global m = m_
    global ϕ = ϕ_
    global xds = xds_
    nothing
end


function calculateclicked()
    ms = speakermix2mat(app[:ent03][:get]())
    mm = microphonemix2mat(app[:ent04][:get]())
    mm = mm[:,:]
    fs = parse(Int, app[:ent00][:get]())
    n = parse(Int, app[:ent01][:get]())
    atten = parse(Float64, app[:ent02][:get]())

    f0 = parse(Float64, app[:ent3][:get]())
    fx = parse(Float64, app[:ent4][:get]())
    f1 = parse(Float64, app[:ent5][:get]())

    @info "calculate:" ms mm fs n atten f0 fx f1
    h_ = ImpulseResponse.minimalphase_fir_design(f, m, ϕ, xds, fs, f0, f1, fx, n)
    global h = h_
    nothing
end


function validateclicked()
    ms = speakermix2mat(app[:ent03][:get]())
    mm = microphonemix2mat(app[:ent04][:get]())
    mm = mm[:,:]
    fs = parse(Int, app[:ent00][:get]())
    n = parse(Int, app[:ent01][:get]())
    atten = parse(Float64, app[:ent02][:get]())

    f0 = parse(Float64, app[:ent6][:get]())
    f1 = parse(Float64, app[:ent7][:get]())

    @info "validate:" ms mm fs n atten f0 f1
    ImpulseResponse.minimalphase_fir_verification(h, ms, mm, fs, f0, f1, atten, n)
end


function saveclicked()
    @info "save to file" h
    h5write(app[:ent8][:get](), app[:ent9][:get](), convert(Matrix{Float64},h[:,:]))
end


app[:mainloop]()