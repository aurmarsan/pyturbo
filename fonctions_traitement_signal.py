import scipy, numpy

#__________________________________________________________________________________________
def recaler_signal(Signal, SignalRef):
    """resample si les longueurs ne sont pas les memes
    C'est alors mieux quand SignalRef est periodique"""
    if SignalRef.size != Signal.size:
        SignalRef = scipy.signal.resample(SignalRef, Signal.size)
    iroll = numpy.argmax(
        numpy.correlate(numpy.r_[Signal, Signal], SignalRef, mode='valid'))
    return numpy.roll(Signal, -iroll), -iroll
#__________________________________________________________________________________________
