function isoct = isOctave()
% Returns true if function is called from Octave (as opposed to Matlab)

    isoct = exist('OCTAVE_VERSION', 'builtin') ~= 0;

return
