#include "Synth.h"
#include <cmath>
#include <cstddef>
#include <sndfile.h>

Synth::Synth(int samplerate, Envelope env)
    : samplerate(samplerate),
      env(env) {
  envToSamples(env);
}

void Synth::addNote(Note note, std::size_t start, std::size_t end) {
  std::vector<double> sine = genSineWave(note, start, end);

  std::size_t required_size = start + sine.size();

  if (required_size > audio.size()) {
    audio.resize(required_size, 0.0);
  }

  for (std::size_t i = 0; i < sine.size(); i++) {
    audio[start + i] += sine[i];
  }
}
void Synth::create_wav(const char *file_name) {
  // setting up file info
  SF_INFO sfinfo = {};
  sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
  sfinfo.samplerate = this->samplerate;
  sfinfo.channels = 1;

  // normalise audio before writing
  normaliseAudio();

  // create and write to file
  SNDFILE *outfile = sf_open(file_name, SFM_WRITE, &sfinfo);
  sf_write_double(outfile, this->audio.data(), this->audio.size());
  sf_close(outfile);
}

void Synth::normaliseAudio() {
  double max_amp = 0.0;
  for (double s : audio)
    max_amp = std::max(max_amp, std::abs(s));

  if (max_amp > 1.0) {
    for (double &s : audio)
      s /= max_amp;
  }
}

double Synth::getFrequency(Note note) {
  int n = static_cast<int>(note.note) + (note.octave - 4) * 12;
  return 440.0 * pow(2.0, (n - 9) / 12.0);
}

std::vector<double> Synth::genSineWave(Note note, std::size_t start,
                                       std::size_t end) {
  std::size_t duration = end - start;
  std::size_t total_length = duration + release_smp;
  std::vector<double> wave(total_length, 0.0);

  double freq = getFrequency(note);
  const double phase_incr = (2 * M_PI * freq) / samplerate;

  for (std::size_t i = 0; i < total_length; i++) {
    double amplitude = 0.0;

    if (i < attack_smp) {
      // Stage 1: Attack (0.0 to 1.0)
      amplitude = static_cast<double>(i) / attack_smp;
    } else if (i < attack_smp + decay_smp) {
      // Stage 2: Decay (1.0 down to sustain_level)
      double decay_progress = static_cast<double>(i - attack_smp) / decay_smp;
      amplitude = 1.0 - (decay_progress * (1.0 - env.sustain_level));
    } else if (i < duration) {
      // Stage 3: Sustain (steady at sustain_level)
      amplitude = env.sustain_level;
    } else {
      // Stage 4: Release (sustain_level down to 0.0)
      double release_progress = static_cast<double>(i - duration) / release_smp;
      amplitude = env.sustain_level * (1.0 - release_progress);
    }

    wave[i] = sin(phase_incr * i) * amplitude;
  }
  return wave;
}

void Synth::envToSamples(Envelope env) {
  attack_smp = msToSamples(env.attack_ms);
  decay_smp = msToSamples(env.decay_ms);
  release_smp = msToSamples(env.release_ms);
}

int Synth::msToSamples(double ms) {
  return static_cast<int>((ms / 1000.0) * samplerate);
}
