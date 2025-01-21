**Loss of pair formation predates the evolution of male-less society in
*Glyptotermes* termites**

**Nobuaki Mizumoto<sup>1\*</sup>, Simon Hellemans<sup>2</sup>**

1\. Department of Entomology & Plant Pathology, Auburn University,
Auburn, AL, 36849, USA

2\. Okinawa Institute of Science & Technology Graduate University,
Onna-son, Okinawa, 904-0495 Japan

\*: Correspondence: <nzm0095@auburn.edu>

**Abstract**

Parthenogenesis and the loss of males have occurred repeatedly across
diverse organisms. Asexually-reproducing lineages are not usually
associated with social animals that exhibit biparental care because such
care is inherently linked to the behavioral sequence of mate pairing and
sexual reproduction. The male-less lineages of the termite,
*Glyptotermes nakajimai*, provides a rare opportunity to study how
sexual reproduction can be lost in social animals with parental care.
Here we demonstrate that modification of the mate-pairing process
predated the evolution of asexual lineages. Termite colonies are
typically founded by a mating pair, with many species forming a tandem
courtship pair while searching for a nest site. Our comparative analysis
of tandem running in *Glyptotermes* termites revealed that two related
species, *G. fuscus* and *G. satsumensis*, exhibited both female-leader
and male-leader tandem runs, estimated to be the ancestral state in this
genus. On the other hand, tandem running was rare and ephemeral in both
sexual and asexual lineages of *G. nakajimai*. These results suggest
that *G. nakajimai* employs an alternative colony foundation strategy,
as further supported by their colony structures. Our study highlights
the tight link between the evolution of asexuality and behavioral
preadaptation, contributing to our understanding of the evolution of
complex phenotypes.

**Keywords:** Asexual reproduction, Movement coordination, Parental
care, Same-sex sexual behavior, Social insects

**Introduction**

Sexual reproduction is norm for multicellular organisms in spite of two
fold costs. Many studies have revealed the advantages of sexual
reporoduction in short and long period. Furthermore, in established
sexual species, there is constrains to prevent from the evolution of
asexual by maintaining low transration rate from sexual to asexual
reproduction \[1\]. Therefore, although asexual lineages evolved across
diversity of taxa, the evolution of asexual population is relatively
rare. In addition to genetic mechanisms, specific behaviors of sexual
organisms can prevent the evolution of asexual linegaes, such as sexual
conflicts (\[2\] but see \[3\]). In social animals with parental care,
the evolution of asexuality is challenging because sexual reproduction
is associated with mate pairing and biparental care, where males
contributes to fitness more than sperm \[4\].

Termites evolved from the ancestor of subsocial lineages, and thus also
colonies start from subsocial pairs. During colony foundation processes,
males significantly contribute to the offspring both energetically and
through physical labor (refs), and thus nest establishment as a single
termite is not usually successful (refs). In this sense, even if
termites have parthenogenetic ability, they need a partner for colony
foundation. A species of drywood termite, *Glyptotermes nakajimai*
Morimoto (Isoptera: Kalotermitidae) \[5\], provides a unique opportunity
to study the evolution of asexual lineages in animals with biparental
care. In this species, all colonies are comprised only of females (i.e.,
all-female asexual societies) in several popualtions \[6\].

There are two different potential behavioral preadaptations that enable
the evolution of a male-less colony foundation in termites. First,
colony establishment by female-female pairs after same-sex tandem runs.
For example, in *Reticulitermes* termites, same-sex tandem runs function
as heterosexual tandem runs \[7\], and female-female pairs start nest
with parthenogenesis \[8–10\], although they cannot grow to the mature
colony (ref). If the ancestor of G. nakajimai has strong tandem running
behavior with same-sex pairing, that can facilitate the evolution of
male-less colony foundation. Second,

pleometrosis

1\. Same-sex tandem runs, female-female tandem runs

2\. Colony foundation by multiple individuals, pairing not based on
tandem running.

**Memo for methods**

*Termite collection*

We collected all termite colonies with a piece of nesting wood from the
filed. We collected three colonies of *G. fuscus* (one in January 2021
and one in March 2022 in Nago, Okinawa; one in March 2023, Iriomote Is.,
Okinawa), three colonies of *G. satsumensis* in March 2021 (two in
Minamiosumi, Kagoshima, one in Kushima, Miyazaki), and four colonies of
*G. nakajimai* (two in March 2021, Wakasa, Fukui, one in April 2023,
Tokunoshima Is. Kagoshima, one in March 2021 in Cape Toi, Miyazaki). For
*G. nakajimai*, samples from Fukui and Tokunoshima Is. were sexual
\[11\], while the sample from Cape Toi was asexual \[6\]. The field
collection was performed before the swarming seasion; each colony
contained nymphs but not alates. All colonies were maintained within the
nesting wood at 22°C until the experiments. Before each experiment, we
transferred nests to a room at 27 °C, which promoted alates to emerge
and fly. Alates were then collected and separated individually. Tandem
running behavior happens after termites shed their wings. We used
individuals that shed their wings by themselves within 12 h.

Alates were then collected, separated by sex, and color-marked with one
dot of paint (PX-20; Mitsubishi) on the abdomen to distinguish sex
identities.

*Behavioral observations*

We introduced a female-male pair of termite dealates (female-female pair
for asexual population of *G. nakajimai*) to the experimental arena,
consisting of a petri dish (φ = 90 mm) covered with a layer of moistened
plaster. All pairs were prepared using nest mates. We recorded their
behavior up to 60 minutes at the rate of 30 frames per second. All the
videos were cropped to 1200x1200 pix to only include the arena in the
frame before the video analysis. In total, we observed 21 pairs of *G.
stsumensis* (340:16, 347:2, JP21-06:3), 46 pairs of *G. fuscus* (21A:16,
G05:18, NM2325:12), 25 pairs of *G. nakajimai* sexual populations
(356:6, 367:3, NM2344:16), and 15 pairs of asexual populations.

All videos were analyzed using SLEAP v 1.4.0 \[12\] to estimate the
movement of body parts of each individual. We used a 6-node skeleton:
antenna tips (LR), head (middle of mouth parts), head-pronotum boundary,
body center, abdomen-tip, and a dot of color painted marker. We built a
model for one species and then used it as a starting point to build
another for the next species sequentially. First, we labeled 342
individuals from 23 videos for training in *G. satsumensis*. We trained
a U-Net-based model with a multi-animal top-down approach, with a
receptive field size of 76 pixels for the centroid and 156 pixels for
the centered instance, on Nvidia GeForce RTX 4090, where augmentation is
done by rotating images from -180 to 180 degrees. The mean Average
Precision (mAP) and mean Average Recall (mAR) of this model were 0.36
and 0.49, respectively. While tracking after the inference, we used the
instance similarity method with the greedy matching method. All pose
estimation data were converted to HDF5 files for further analysis. In
*G. fuscus*, TBA.

We used Python to format all HDF5 files for further analysis and
converted them into FEATHER files for analysis in R \[13\]. We employed
a linear interpolation method to address missing values in the dataset.
After scaling all data from pixels to mm (1200 pixels = arena size), we
used a median filter with a kernel size of 5 to reduce noise.

To compare tandem running behaviors among species, we automatically
determined that pairs were in tandem or not based on the postures and
spatial position of partners. First, we regarded two individuals were in
interaction when the distance between body centers of partners was less
than two body lengthes, based on the frequency distribution of this
distance (Fig. 1C). In this process, we ignored the short interaction
events or non-interaction events less than 2 seconds to smooth the data.
Second, during interactions, we classified termite heading orientation
as female-leader and male-leader. We obtained heading directions of
females and males as vector from abdomen tips to head front. Then a pair
was in female-leader when male was behind relative to female heading
direction, and female was front relative to male heading direction, and
vice versa (Fig. 1AB). If a pair spent in female-leader position for
more than half of the time during an interaction event, we regarded that
the interaction event was female-leader tandem runs. This classified all
frames into female-leader tandem, male-leader tandem, other interactions
(including tandem runs where they switch leader-follower roles), and
non-interactions. We obtrained the traveled distance for which the
leader walked during each tandem running events. Then we compared this
traveled distance, using mixed-effects Cox models, with species being
treated as a fixed effect and each pair id as a random effect. We used
coxme() function in the coxme package in R \[14\]. Note that we used
distance instead of duration to evaluate how much tandem running pair
could explore the envronments by removing pausing time during
interactions.

<img src="media/image1.png" style="width:4.62014in;height:3.72222in" />

> **Figure 1.** Spatial positioning between partners in *Glyptotermes*
> termites. (A) Comparison of the relative position of the partner,
> given that female (left) or male (right) heading towards the top at
> the center. Simplified phylogenetic relationship based on \[6\] is
> also provided. (B) Distributions of the partner's position relative to
> the female's heading direction in angles when the pair is within 2
> body lengths. (C) Distributions of the distance between partners.

<img src="media/image2.png" style="width:4.91736in;height:3.56667in" />

> **Figure 2.** Tandem running behavior of each species. (A) Proportion
> of time in each state during observation. Each bar represents one
> pair. (B-C) Interspecific comparison of the traveled distance during
> each tandem running event.

**Idea of another result to be shown.**

I want to show that the colonies of G. nakajimai include many
reproductives (not a pair) \[6\], while colonies of G. fuscus and G.
satsumensis often have just monogamous pairs (often physogastric). This
is consistent with my observation, but my data is minimal as I have not
recorded it properly. The former is shown in the paper \[6\], but we do
not have published information on the latter.

Here are some possible approaches:

Based on this project report (you can read this in English from this
link:
<https://koara.lib.keio.ac.jp/xoonips/modules/xoonips/detail.php?koara_id=2018000005-20180253>),
G. satsumensis keeps monogamous pairing. However, this report is not
published (or looks like will not be published in the future). We may
reach out Hayashi-san to ask if he has any data about colony structures
in G. satsumensis. If so, we can ask if he can join this paper.

We can also reach out Yashiro-san if he have data.

Do you have any thoughts/ideas?

**Idea of Discussion**

Although most theoretical studies of the evolution of sexuality have
actually acknowledged that nonrandom mating or parental care may
infuence the outcome of their models, the importance of these phenomena
has always been minimized \[4\].

**Acknowledgments**

We thank Kensei Kikuchi and Esra Kaymak for field collection, Aoi
Mizumoto for assist in video recording. The work was supported by a JSPS
Research Fellowship for Young Scientists CPD to NM (20J00660), a
Grant-in-Aid for Early-Career Scientists (21K15168) to NM.

HATCH project number.

**References**
