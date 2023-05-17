from loica import *
from flapjack import *
import flapjack as fj
import getpass


user = input()
passwd = getpass.getpass()
flap = fj.flapjack.Flapjack('flapjack.rudge-lab.org:8000')
flap.log_in(username=user, password=passwd)

study = flap.get('study', name='SC loica source rate 100, deg rate 1', description='SC study for demonstrating source rate 100, deg rate 1')
dna = flap.get('dna', name='source')
vector = flap.get('vector', name='source', dnas=dna.id)
gfp = flap.get('signal', name='GFP', color='green', description='fluorescent')
media = flap.get('media', name='loica', description='Simulated loica media')
strain = flap.get('strain', name='loica', description='Loica test strain')

reporter = Reporter(name='GFP', color='green', init_concentration=0, signal_id=gfp.id[0])

source = Source(output=reporter)

source.characterize_stochastic(flap,
    vector=vector.id,
    media=media.id,
    strain=strain.id,
    signal=gfp.id,)
