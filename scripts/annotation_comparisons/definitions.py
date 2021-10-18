
from dataclasses import dataclass


@dataclass(init=False)
class LRC:
    '''
    The 'low-resolution' consensus labels (as simple as it gets)
    '''    
    nucleus = 'nucleus'
    organelles = 'organelles'
    cytoplasm = 'cytoplasm'


@dataclass(init=False)
class HRC:
    '''
    The 'high-resolution' consensus labels
    (to which are mapped the OC, HPA, and yeast annotations)
    '''    
    nuclear_membrane = 'nuclear membrane'
    nucleolus = 'nucleolus'
    nuclear_domains = 'nuclear domains'
    
    # OC-HPA only
    nucleoplasm = 'nucleoplasm'

    cytoplasm = 'cytoplasm'    
    membrane = 'membrane'
    er = 'ER'
    mitochondria = 'mitochondria'
    
    # OC-yeast only
    punctate = 'vesicles and Golgi'

    # OC-HPA only
    vesicles = 'vesicles'
    golgi = 'Golgi'
    centrosome = 'centrosome'
    cytoskeleton = 'cytoskeleton'


class OC:
    nucleoplasm = 'nucleoplasm'
    nucleolus_fc_dfc = 'nucleolus_fc_dfc'
    nucleolus_gc = 'nucleolus_gc'
    nuclear_punctae = 'nuclear_punctae'
    chromatin = 'chromatin'
    nuclear_membrane = 'nuclear_membrane'
    cytoplasmic = 'cytoplasmic'
    cytoskeleton = 'cytoskeleton'
    membrane = 'membrane'
    cell_contact = 'cell_contact'
    focal_adhesions = 'focal_adhesions'
    vesicles = 'vesicles'
    centrosome = 'centrosome'
    mitochondria = 'mitochondria'
    golgi = 'golgi'
    er = 'er'


class HPA:
    nucleoplasm = 'Nucleoplasm'
    nucleoli_rim = 'Nucleoli rim'
    nuclear_speckles = 'Nuclear speckles'
    nuclear_membrane = 'Nuclear membrane'
    nucleoli = 'Nucleoli'
    nucleoli_fibrillar_center = 'Nucleoli fibrillar center'
    nuclear_bodies = 'Nuclear bodies'
    plasma_membrane = 'Plasma membrane'
    cytosol = 'Cytosol'
    mitochondria = 'Mitochondria'
    endoplasmic_reticulum = 'Endoplasmic reticulum'
    golgi_apparatus = 'Golgi apparatus'
    centrosome = 'Centrosome'
    peroxisomes = 'Peroxisomes'
    endosomes = 'Endosomes'
    lysosomes = 'Lysosomes'
    vesicles = 'Vesicles'
    lipid_droplets = 'Lipid droplets'
    cell_junctions = 'Cell Junctions'
    focal_adhesion_sites = 'Focal adhesion sites'
    microtubules = 'Microtubules'
    intermediate_filaments = 'Intermediate filaments'
    actin_filaments = 'Actin filaments'
    centriolar_satellite = 'Centriolar satellite'
    cytoplasmic_bodies = 'Cytoplasmic bodies'


class Yeast:
    nucleus = 'nucleus'
    nucleolus = 'nucleolus'
    nuclear_periphery = 'nuclear periphery'
    cytosol = 'cytosol'
    vacuole_membrane = 'vacuole membrane'
    vacuole = 'vacuole'
    punctate = 'punctate'
    cell_periphery = 'cell periphery'
    mitochondria = 'mitochondria'
    er = 'ER'


low_res_oc_yeast_labels = [
    {
        'consensus': LRC.nucleus,
        'oc': (
            OC.nucleoplasm, 
            OC.nucleolus_fc_dfc, 
            OC.nucleolus_gc, 
            OC.nuclear_punctae, 
            OC.chromatin,
            OC.nuclear_membrane,
        ),
        'yeast': (Yeast.nucleus, Yeast.nucleolus, Yeast.nuclear_periphery),
    },{
        'consensus': LRC.cytoplasm,
        'oc': (
            OC.cytoplasmic, 
            # OC.cytoskeleton,
        ),
        'yeast': (Yeast.cytosol,),
    },{
        'consensus': LRC.organelles,
        'oc': (
            OC.membrane, 
            OC.cell_contact, 
            OC.focal_adhesions,
            OC.vesicles,
            OC.centrosome, 
            OC.golgi, 
            OC.mitochondria, 
            OC.er, 
        ),
        'yeast': (
            Yeast.vacuole_membrane, 
            Yeast.vacuole, 
            Yeast.punctate,
            Yeast.cell_periphery, 
            Yeast.mitochondria, 
            Yeast.er,
        ),
    },
]


high_res_oc_yeast_labels = [
    {
        'consensus': HRC.nuclear_domains,
        'oc': (OC.nucleoplasm, OC.nuclear_punctae, OC.chromatin,),
        'yeast': Yeast.nucleus
    },{
        'consensus': HRC.nucleolus,
        'oc': (OC.nucleolus_fc_dfc, OC.nucleolus_gc,),
        'yeast': Yeast.nucleolus
    },{
        'consensus': HRC.nuclear_membrane,
        'oc': OC.nuclear_membrane,
        'yeast': Yeast.nuclear_periphery,
    },{
        # nb this category does not exist in the OC-HPA mapping
        'consensus': HRC.punctate,
        'oc': (OC.vesicles, OC.golgi, OC.centrosome),
        'yeast': (Yeast.vacuole_membrane, Yeast.vacuole, Yeast.punctate)
    },{
        'consensus': HRC.cytoplasm,
        'oc': OC.cytoplasmic,
        'yeast': Yeast.cytosol,
    },{
        'consensus': HRC.membrane,
        'oc': (OC.membrane, OC.cell_contact, OC.focal_adhesions),
        'yeast': Yeast.cell_periphery,
    },{
        'consensus': HRC.mitochondria,
        'oc': OC.mitochondria,
        'yeast': Yeast.mitochondria,
    },{
        'consensus': HRC.er,
        'oc': OC.er,
        'yeast': Yeast.er,
    },
]

low_res_oc_hpa_labels = [
    {
        'consensus': LRC.nucleus,
        'oc': (
            OC.nucleoplasm, 
            OC.nucleolus_fc_dfc, 
            OC.nucleolus_gc, 
            OC.nuclear_punctae, 
            OC.chromatin,
            OC.nuclear_membrane,
        ),
        'hpa': (
            HPA.nucleoplasm, 
            HPA.nucleoli, 
            HPA.nucleoli_fibrillar_center, 
            HPA.nucleoli_rim, 
            HPA.nuclear_speckles, 
            HPA.nuclear_bodies,
            HPA.nuclear_membrane,
        ),
    },{
        # not clear whether 'cytoskeleton' belongs in this category (or any low-res category)
        'consensus': LRC.cytoplasm,
        'oc': (
            OC.cytoplasmic, 
            # OC.cytoskeleton,
        ),
        'hpa': (
            HPA.cytosol,
            # HPA.microtubules, HPA.actin_filaments, HPA.intermediate_filaments,
        ),
    },{
        'consensus': LRC.organelles,
        'oc': (
            OC.membrane, 
            OC.cell_contact, 
            OC.focal_adhesions,
            OC.mitochondria, 
            OC.golgi, 
            OC.er, 
            OC.centrosome, 
            OC.vesicles
        ),
        'hpa': (
            HPA.plasma_membrane, 
            HPA.cell_junctions, 
            HPA.focal_adhesion_sites, 
            HPA.mitochondria, 
            HPA.golgi_apparatus, 
            HPA.endoplasmic_reticulum, 
            HPA.centrosome, 
            HPA.lipid_droplets, 
            HPA.peroxisomes, 
            HPA.lysosomes, 
            HPA.endosomes, 
            HPA.vesicles
        ),
    },
]


high_res_oc_hpa_labels = [
    {
        'consensus': HRC.nuclear_domains,
        'oc': (OC.nuclear_punctae, OC.chromatin),
        'hpa': (HPA.nuclear_speckles, HPA.nuclear_bodies),
    },{
        'consensus': HRC.nucleolus,
        'oc': (OC.nucleolus_fc_dfc, OC.nucleolus_gc,),
        'hpa': (HPA.nucleoli, HPA.nucleoli_fibrillar_center, HPA.nucleoli_rim),
    },{
        'consensus': HRC.nucleoplasm,
        'oc': OC.nucleoplasm,
        'hpa': HPA.nucleoplasm,
    },{
        'consensus': HRC.nuclear_membrane,
        'oc': OC.nuclear_membrane,
        'hpa': HPA.nuclear_membrane,
    },{
        'consensus': HRC.cytoplasm,
        'oc': OC.cytoplasmic,
        'hpa': HPA.cytosol,
    },{
        'consensus': HRC.mitochondria,
        'oc': OC.mitochondria,
        'hpa': HPA.mitochondria,
    },{
        'consensus': HRC.golgi,
        'oc': OC.golgi,
        'hpa': HPA.golgi_apparatus,
    },{
        'consensus': HRC.er,
        'oc': OC.er,
        'hpa': HPA.endoplasmic_reticulum,
    },{
        'consensus': HRC.centrosome,
        'oc': OC.centrosome,
        'hpa': HPA.centrosome,
    },{
        'consensus': HRC.membrane,
        'oc': (OC.membrane, OC.cell_contact, OC.focal_adhesions),
        'hpa': (HPA.plasma_membrane, HPA.cell_junctions, HPA.focal_adhesion_sites),
    },{
        'consensus': HRC.cytoskeleton,
        'oc': OC.cytoskeleton,
        'hpa': (HPA.microtubules, HPA.actin_filaments, HPA.intermediate_filaments),
    },{
        'consensus': HRC.vesicles,
        'oc': OC.vesicles,
        'hpa': (HPA.lipid_droplets, HPA.peroxisomes, HPA.lysosomes, HPA.endosomes, HPA.vesicles),
    },
]

# the order of the consensus labels in the sankey diagram
# (this is manually defined for all single labels and the most common label sets,
# for both the yeast and HPA consensus labels)
label_orders = {
    'high': [
        (HRC.nuclear_membrane,),
        (HRC.nucleolus,),
        (HRC.nuclear_domains,),
        (HRC.nucleoplasm, HRC.nuclear_domains,),
        (HRC.nucleoplasm,),
        (HRC.nucleolus, HRC.punctate),

        (HRC.cytoplasm, HRC.nucleoplasm),
        (HRC.cytoplasm, HRC.nuclear_domains,),
        (HRC.cytoplasm, HRC.nucleolus),
        
        (HRC.cytoplasm,),
        
        (HRC.cytoplasm, HRC.vesicles),
        (HRC.vesicles,),
        
        (HRC.vesicles, HRC.er),
        (HRC.er,),

        (HRC.vesicles, HRC.membrane),
        (HRC.cytoplasm, HRC.membrane),
        (HRC.membrane,),

        (HRC.punctate,),
        (HRC.cytoplasm, HRC.punctate),

        (HRC.cytoplasm, HRC.nuclear_domains, HRC.punctate),
        (HRC.cytoplasm, HRC.membrane, HRC.punctate),
        
        (HRC.punctate, HRC.er),
        (HRC.vesicles, HRC.golgi),

        (HRC.golgi,), 
        (HRC.mitochondria,), 
        (HRC.centrosome,),
        (HRC.cytoskeleton,),
    ],

    'low': [
        (LRC.nucleus,),
        (LRC.nucleus, LRC.cytoplasm),
        (LRC.cytoplasm,),
        (LRC.organelles, LRC.cytoplasm),
        (LRC.organelles,),
        (LRC.nucleus, LRC.organelles),
        (LRC.nucleus, LRC.organelles, LRC.cytoplasm)
    ]
}


@dataclass(init=False)
class ColorsFig7B:
    '''
    Manually-defined colors from the localization hierarchy in Fig 7B of the preprint
    '''    
    # nuclear categories
    orange = '#f89f3d'
    brown = '#714d2d'
    dbrown = '#404041'
    yellow = '#fee044'
    pink = '#f7aab3'
    lbrown = '#b18961'
    red = '#cd5a53'

    # organelles
    purple = '#8c53a2'
    lbluepurple = '#a5c1e5'
    dbluepurple = '#2a3890'
    lblue = '#45b3e4'
    bluegreen = '#48969b'
    dblue = '#3d86bd'

    # cytoplasm
    graygreen = '#a6b7a2'
    lime = '#dee78d'
    dgreen = '#00793e'
    ddgreen = '#254c32'
    green = '#98cb57'

# colors for each consensus category that match Fig 7B of the preprint
consensus_label_colors_fig7b = {
    'low': {
        LRC.nucleus: ColorsFig7B.brown,
        LRC.organelles: ColorsFig7B.lblue,
        LRC.cytoplasm: ColorsFig7B.green,
    },
    'high': {
        HRC.nuclear_membrane: ColorsFig7B.yellow,
        HRC.nucleolus: ColorsFig7B.brown,
        HRC.nuclear_domains: ColorsFig7B.lbrown,
        HRC.nucleoplasm: ColorsFig7B.dbrown,
        HRC.cytoplasm: ColorsFig7B.green,
        HRC.punctate: ColorsFig7B.lblue,
        HRC.membrane: ColorsFig7B.dblue,
        HRC.vesicles: ColorsFig7B.lblue,
        HRC.golgi: ColorsFig7B.dbluepurple,
        HRC.mitochondria: ColorsFig7B.lbluepurple,
        HRC.centrosome: ColorsFig7B.bluegreen,
        HRC.cytoskeleton: ColorsFig7B.bluegreen,
        HRC.er: ColorsFig7B.purple,
    }
}


@dataclass(init=False)
class ColorsHSV:
    '''
    Manually-defined colors from the hard-coded HSV-like colormap
    This version is hard-coded for reproducibility 
    between sankey diagrams with different subsets of label sets
    '''    
    red = '#f44336'
    pink = '#e91e63'
    lpurple = '#9c27b0'
    dpurple = '#673ab7'
    dblue = '#3f51b5'
    blue = '#03a9f4'
    lblue = '#00bcd4'
    dgreen = '#009688'
    green = '#4caf50'
    lgreen = '#8bc34a'
    llgreen = '#cddc39'
    yellow = '#ffeb3b'
    orange = '#ffc107'

consensus_label_colors = {
    'low': {
        LRC.nucleus: ColorsHSV.lgreen,
        LRC.organelles: ColorsHSV.red,
        LRC.cytoplasm: ColorsHSV.blue,
    },
    'high': {
        HRC.nuclear_membrane: ColorsHSV.yellow,
        HRC.nucleolus: ColorsHSV.llgreen,
        HRC.nuclear_domains: ColorsHSV.lgreen,
        HRC.nucleoplasm: ColorsHSV.dgreen,
        HRC.cytoplasm: ColorsHSV.lblue,
        HRC.punctate: ColorsHSV.dblue,
        HRC.vesicles: ColorsHSV.dpurple,
        HRC.membrane: ColorsHSV.lpurple,
        HRC.er: ColorsHSV.pink,
        HRC.golgi: ColorsHSV.red,
        HRC.mitochondria: ColorsHSV.red,
        HRC.centrosome: ColorsHSV.orange,
        HRC.cytoskeleton: ColorsHSV.orange,
    }
}