#include <cmath>
#include <iostream>
#include <sndfile.h>
#include <vector>

enum class Notes { C, Cs, D, Ds, E, F, Fs, G, Gs, A, As, B };

struct Note {
  Notes note;
  int octave;
};

enum class Mode { Major, Minor };

struct Scale {
  Notes root;
  Mode mode;
};

class Synth {
public:
  Scale scale;
  std::vector<double> audio;
  int samplerate;

  Synth(Scale scale, int samplerate) {
    this->scale = scale;
    this->samplerate = samplerate;
  }

  void add_note(struct Note note, std::size_t start, std::size_t end) {
    // if end sample is beyond current audio size, resize audio
    if (end >= audio.size()) {
      this->audio.resize(end, 0.0);
    }
    std::vector<double> sine = gen_sine_wave(note, end - start);
    // add beep to audio at specified position
    for (std::size_t i = start; i < end; i++) {
      this->audio[i] += sine[i - start];
    }
  }

  void create_wav() {
    // setting up file info
    SF_INFO sfinfo = {};
    sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
    sfinfo.samplerate = this->samplerate;
    sfinfo.channels = 1;

    // normalise audio before writing
    normalise_audio();

    // create and write to file
    SNDFILE *outfile = sf_open("output.wav", SFM_WRITE, &sfinfo);
    sf_write_double(outfile, this->audio.data(), this->audio.size());
    sf_close(outfile);
  }

private:
  void normalise_audio() {
    double max_amp = 0.0;
    for (double s : audio)
      max_amp = std::max(max_amp, std::abs(s));

    if (max_amp > 1.0) {
      for (double &s : audio)
        s /= max_amp;
    }
  }

  double get_frequency(Note note) {
    int n = static_cast<int>(note.note) + (note.octave - 4) * 12;
    return 440.0 * pow(2.0, (n - 9) / 12.0);
  }

  std::vector<double> gen_sine_wave(Note note, std::size_t num_samples) {
    // initialise beep vector
    std::vector<double> wave(num_samples, 0.0);
    double freq = get_frequency(note);
    // find phase increment
    const double phase_incr = (2 * M_PI * freq) / this->samplerate;
    // fill beep vector with sine wave samples
    for (std::size_t i = 0; i < wave.size(); i++) {
      wave[i] = sin(phase_incr * i);
    }
    return wave;
  }
};

int main() {
  SF_INFO sfinfo = {};
  SNDFILE *sf = sf_open("../assets/ball_bounce.wav", SFM_READ, &sfinfo);
  if (!sf) {
    std::cerr << "Error: " << sf_strerror(NULL);
    throw std::runtime_error("Failed to open audio file");
  }

  // limit file size to ~30 minutes at 44.1kHz
  if (sfinfo.frames / sfinfo.channels > 158760000) {
    std::cerr << "file too large, exiting";
    throw std::runtime_error("File too large");
  }

  // read audio to vector
  std::vector<double> samples(sfinfo.frames * sfinfo.channels);
  sf_readf_double(sf, samples.data(), sfinfo.frames);
  sf_close(sf);

  // convert to mono signal
  std::vector<double> mono;
  if (sfinfo.channels == 2) {
    mono.resize(sfinfo.frames);
    for (sf_count_t i = 0; i < sfinfo.frames; i++) {
      mono[i] = 0.5f * (samples[i * 2] + samples[i * 2 + 1]);
    }
  } else {
    mono = samples;
  }

  // create synth
  Scale scale = {Notes::C, Mode::Major};
  Synth synth(scale, sfinfo.samplerate);
  // c major chord

  synth.add_note({Notes::C, 4}, 0, 100000);
  synth.add_note({Notes::E, 4}, 0, 100000);
  synth.add_note({Notes::G, 4}, 0, 100000);
  synth.add_note({Notes::B, 4}, 0, 100000);
  synth.create_wav();

  return 0;
}
