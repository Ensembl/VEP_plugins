=head1 LICENSE

 Copyright (c) 1999-2012 The European Bioinformatics Institute and                                                   
 Genome Research Limited.  All rights reserved.                                                                      

 This software is distributed under a modified Apache license.                                                       
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html                                                               

=head1 CONTACT                                                                                                       

 Will McLaren <wm2@ebi.ac.uk>
    
=cut

=head1 NAME

 Draw

=head1 SYNOPSIS

 mv VEPbrowse.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin VEPbrowse --individual JohnDoe

=head1 DESCRIPTION

 A VEP plugin that generates HTML with an overview of an individual's mutations
 across the genome, as well as a transcript-specific sequence view. This plugin
 is under development, please expect bugs!!!

=cut

package VEPbrowse;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(MAX_DISTANCE_FROM_TRANSCRIPT overlap);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(align_seqs);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::SeqUtils;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

my $ktbs;

sub new {
    my $class = shift;

    my $self = $class->SUPER::new(@_);
    
    # configure
    my @params = @{$self->params};
    
    # initialize cache
    $self->{cache} = {};
    $self->{has_cache} = 1;
    
    # force various options
    $self->{config}->{prefetch} = 1;
    $self->{config}->{hgnc} = 1;

    return $self;
}

sub version {
    return '2.5';
}

sub feature_types {
    return ['Transcript'];
}

sub variant_feature_types {
    return ['VariationFeature'];
}

sub get_header_info {
    return {};
}

# the destroy method will be run as the script closes
# so we can use this to build the summary HTML page
sub DESTROY {
    my $self = shift;
    
    open OUT, ">index.html";
    
    if(defined($self->{karyotype_bands})) {
        my ($canvas_width, $canvas_height) = (1000, 500);
        
        my $html =<<END;
<!DOCTYPE html>
<html>
<head>
</head>
<body onload="init();">

<div id="canvas_holder" style="border:1px solid #c3c3c3; width:${canvas_width}px;">
    <canvas id="karyotype" width="$canvas_width" height="$canvas_height">
    Your browser does not support the HTML5 canvas element.
    </canvas>
</div>
<div id="info" style="width:${canvas_width}px;height:200px;overflow:auto;border:1px solid #c3c3c3; padding: 5px; font-family: arial, sans-serif;">
</div>

<script type="text/javascript">

var c, ctx, tool, genes, info;

var select_from, select_to;
var canvas_width, canvas_height;
var x_scale, x_off;

function init() {
    c = document.getElementById("karyotype");
    ctx = c.getContext("2d");
    
    info = document.getElementById("info");
    
    canvas_width  = $canvas_width;
    canvas_height = $canvas_height;
    
    tool = new tool_pencil();
    
    drawKaryotype();
    
    c.addEventListener('mousedown', ev_canvas, false);
    c.addEventListener('mousemove', ev_canvas, false);
    c.addEventListener('mouseup',   ev_canvas, false);
}
END
        
        $html .= $self->karyotype_jscript();
        $html .= $self->render_karyotype_html($canvas_width, $canvas_height);
        $html .= "\n</script>\n</body>\n</html>";

        print OUT $html;
    }
    
    else {
        print OUT "<ul>";
        foreach my $gene(keys %{$self->{genes}}) {
            print OUT "<li>$gene<ul>\n";
            foreach my $tr(@{$self->{genes}->{$gene}->{transcripts}}) {
                print OUT "<li><a href='".$tr->{transcript}.".html'>$tr->{transcript}</a>";
                print OUT " Mismatches: ".$tr->{mismatches}." Length ratio: ".$tr->{length_ratio};
                print OUT "</li>";
            }
            print OUT "</ul></li>";
        }
        print OUT "</ul>";
    }
    
    close OUT;
}

sub run {
    my ($self, $tva) = @_;
    
    my $tr = $tva->transcript;
    my $vf = $tva->base_variation_feature;
    
    my $tr_stable_id = $tr->stable_id;
    my $individual = $vf->{individual};
    
    $self->get_karyotype_bands($self->config->{reg}->get_adaptor($self->config->{species}, 'core', 'slice')->db);
    
    return {} unless defined($vf->{last_in_transcript}->{$tr_stable_id});
    return {} unless defined $tr->{vfs} && scalar @{$tr->{vfs}};
    
    my $exons = $tr->get_all_Exons();
    return {} unless defined $exons && scalar @$exons;
    $_->{slice} ||= $tr->{slice} for @$exons;
    
    my $mapper = $tr->{_variation_effect_feature_cache}->{mapper};
    return {} unless defined $mapper;
    
    my $codon_table = $tr->{_variation_effect_feature_cache}->{codon_table};
    
    return {} unless defined $tr->translation;
    
    my $translateable_seq = $tr->{_variation_effect_feature_cache}->{translateable_seq};
    
    # check if this transcript has any variants that change its sequence
    my $ok = 0;
    
    my $worst_con;
    
    foreach my $vf(@{$tr->{vfs}}) {
        foreach my $tv(@{$vf->get_all_TranscriptVariations([$tr])}) {
            foreach my $tva(@{$tv->get_all_TranscriptVariationAlleles}) {
                foreach my $oc(@{$tva->get_all_OverlapConsequences}) {
                    my $so_term = $oc->SO_term;
                    $worst_con = $oc if !defined($worst_con) || $oc->rank < $worst_con->rank;
                    
                    if(
                        $so_term eq 'missense_variant' ||
                        $so_term eq 'stop_gained' ||
                        $so_term eq 'stop_lost' ||
                        $so_term eq 'transcript_ablation' ||
                        $so_term eq 'frameshift_variant' ||
                        $so_term eq 'initiator_codon_variant' ||
                        $so_term =~ /inframe_codon/
                    ) {
                        $ok = 1;
                        last;
                    }
                }
            }
        }
    }
    
    return {} unless $ok;
    
    #$DBsingle = 1;
    
    #print STDERR "\nDoing $tr_stable_id $individual\n";
    #return {} unless $tr_stable_id eq 'ENST00000343703';
    
    #print STDERR "Getting ref features\n";
    
    # get ref features
    my @ref_features = @{$self->get_ref_transcript_features($tr)};
    
    # check ref features
    die("ERROR: Ref features borked\n") unless grep {$_->{t} eq 'coding' && $_->{seq}} @ref_features;
    
    #print STDERR "Getting ind features\n";
    
    # get aligned individual features
    # this also re-aligns the reference
    my $features_arrayref = $self->get_individual_features($tr, \@ref_features);
    
    my $var_features = pop @$features_arrayref;
    
    #print STDERR "Getting translations\n";
    
    # get peptides
    my @peptides = map {$self->translate_features($_, $codon_table)} (\@ref_features, @$features_arrayref);
    
    
    #print STDERR "Aligning ref peptide\n";
    
    # align peptides to features
    $self->align_pep_to_features($peptides[0], \@ref_features);
    
    
    #print STDERR "Aligning ind peptides\n";
    my $i = 1;
    $self->align_pep_to_features($peptides[$i++], $_) for @$features_arrayref;
    
    
    #print STDERR "Marking up codons\n";
    
    # mark up codons
    $self->markup_codons($_) for (\@ref_features, @$features_arrayref);
    
    
    #print STDERR "Compiling sequences\n";
    
    # construct seq from features
    my $ref_tr_seq = "";
    $ref_tr_seq .= $_->{marked_seq} for @ref_features;
    
    # compile sequences, ref first, tr then pep
    my @seqs = ($ref_tr_seq);
    
    my $seq = "";
    $seq .= $_->{pep} for @ref_features;
    push @seqs, $seq;
    
    for my $phase(@$features_arrayref) {
        my $seq = "";
        $seq .= $_->{marked_seq} for @{$phase};
        push @seqs, $seq;
        
        $seq = "";
        $seq .= $_->{pep} for @{$phase};
        
        push @seqs, $seq;
    }
    
    $seq = "";
    $seq .= $_->{seq} for @$var_features;
    push @seqs, $seq;
    
    my $seqs = join ",", map {"'$_'"} @seqs;
    
    # compare peptides
    my $ref_pep = $seqs[1];
    
    my ($worst_length_ratio, $worst_mismatches) = (100, 0);
    
    for(my $i=3; $i<=$#seqs; $i+=2) {
        my $pep = $seqs[$i];
        
        # seq_differences
        my ($tp1, $tp2) = ($ref_pep, $pep);
        $tp1 =~ s/\-|\_//g;
        $tp2 =~ s/\-|\_//g;
        
        my $length_ratio = $tp1 ? 100 * (length($tp2) / length($tp1)) : 0;
        $worst_length_ratio = $length_ratio if abs($length_ratio - 100) > abs($worst_length_ratio - 100);
        my $mismatches = 0;
        
        for(my $j=0; $j<length($tp1) && $j<length($tp2); $j+=3) {
            $mismatches++ unless substr($tp1, $j, 3) eq substr($tp2, $j, 3);
        }
        
        $worst_mismatches = $mismatches if $mismatches > $worst_mismatches;
    }
    
    #my @seqs = ([map {$_->{marked_seq}} @ref_features],[map {$_->{pep_seq}} @ref_features]);
    #
    #foreach my $phase(@$features_arrayref) {
    #    push @seqs, [map {$_->{marked_seq}} @{$phase}];
    #    push @seqs, [map {$_->{pep_seq}} @{$phase}];
    #}
    #
    #my $seqs = '[';
    #
    #foreach my $seq(@seqs) {
    #    $seqs .= '['.(join ",", map {"'$_'"} @$seq).'],';
    #}
    #
    #$seqs =~ s/\,$//g;
    #$seqs .= ']';
    
    # log gaps in align
    my @gaps = ('SHIFTMEOFF');
    for my $i(0..(length($ref_tr_seq) - 1)) {
        push @gaps, $i + 1 if substr($ref_tr_seq, $i, 1) eq '-';
    }
    my $gaps = join ",", map {"'$_'"} @gaps;
    
    # parse VFs
    my $vfs = $self->format_vfs($tr);
    
    my $canvas_width  = 1000;
    my $canvas_height = 100;
    
    my $x_scale  = ($canvas_width - 20)  / ($tr->end - $tr->start + 1);
    my $y_scale  = ($canvas_height - 20) / 100;
    my $x_off    = 10;
    
    my $css = get_css();
    
    my $html =<<END;
<!DOCTYPE html>
<html>
<head>
$css
</head>
<body onload="init();">

<div id="header" style="width:${canvas_width}px" class="header">
    <div style="display: inline; float: right;"><img src="http://www.ensembl.org/img/vep_logo.png"></div>
    <h1 style="display: inline; vertical-align: top;">$tr_stable_id</h1><br/>
    <h2 style="display: inline; vertical-align: top;">$individual</h2>
    <div style="font-size:small"><a href="index.html">Back to index</a></div>
</div>

<div id="canvas_holder" style="border:1px solid #c3c3c3; width:${canvas_width}px;">
    <canvas id="transcript" width="$canvas_width" height="$canvas_height">
    Your browser does not support the HTML5 canvas element.
    </canvas>
</div>

<div id="control" style="width:${canvas_width}px" class="control">
    <b>Selection: </b>
    <a href="#" onClick="select_all();">Select all</a>
    |
    <a href="#" onClick="clear_selection();">Clear selection</a>
    |
    <b>Zoom: </b><a href="#" onClick="zoom(0.75);">+</a> / <a href="#" onClick="zoom(1.33);">-</a>
    |
    <b>Move: </b><a href="#" onClick="move(-0.5);">&lt;&lt;</a> / <a href="#" onClick="move(-0.1);">&lt;</a> / <a href="#" onClick="move(0.1);">&gt;</a> / <a href="#" onClick="move(0.5);">&gt;&gt;</a>
    
    <hr class="control_spacer">
    
    <b>Options: </b>
    <input type="checkbox" name="coding_only_checkbox" id="coding_only_checkbox" onclick="coding_only=toggle_value(coding_only);render_seq();" checked=1/>
    <span onClick="var c = document.getElementById('coding_only_checkbox'); c.checked =(! c.checked); coding_only=toggle_value(coding_only);render_seq();">Show coding sequence only</span>
    
    <hr class="control_spacer">
    
    <b>Sequences: </b>
    Reference
    <select name="ref_seqs_select" id="ref_seqs_select" onchange="ref_seqs=this.value;render_seq();">
        <option value="tr_pep">Both</option>
        <option value="tr">Transcript</option>
        <option value="pep">Peptide</option>
        <option value="none">None</option>
    </select>
    |
    Phase 1
    <select name="p1_seqs_select" id="p1_seqs_select" onchange="p1_seqs=this.value;render_seq();">
        <option value="tr_pep">Both</option>
        <option value="tr">Transcript</option>
        <option value="pep">Peptide</option>
        <option value="none">None</option>
    </select>
    |
    Phase 2
    <select name="p2_seqs_select" id="p2_seqs_select" onchange="p2_seqs=this.value;render_seq();">
        <option value="tr_pep">Both</option>
        <option value="tr">Transcript</option>
        <option value="pep">Peptide</option>
        <option value="none">None</option>
    </select>
    
</div>

<div id="sequence" style="width:${canvas_width}px;height:400px;overflow:auto;border:1px solid #c3c3c3">
    Click and drag to select a region of the transcript and show sequence.
</div>
<div id="debug" style="width:${canvas_width}px;height:150px;overflow:auto;border:1px solid #c3c3c3; font-size: small;"><h3>Debug</h3></div>

<script type="text/javascript">

var c, ctx, tool, s, debug_div, canvas_scale;

var select_from, select_to;

var seqs, vfs;
var coding_only, ref_seqs, p1_seqs, p2_seqs;
var seq_div_html;
var canvas_width, canvas_height;
var x_scale, x_off;

var gaps;

function init() {
    c = document.getElementById("transcript");
    ctx = c.getContext("2d");
    
    s = document.getElementById("sequence");
    debug_div = document.getElementById("debug");
    
    canvas_width  = $canvas_width;
    canvas_height = $canvas_height;
    x_scale       = $x_scale;
    x_off         = $x_off;
    
    gaps = new Array($gaps);
    seqs = new Array($seqs);
    vfs = $vfs;
    
    // shift off first element
    // bugfix for strange javascript behaviour when you try to create an array
    // with one numeric element
    gaps.shift();

    drawTranscript();
    
    canvas_scale = 1;
    
    coding_only = 1;
    ref_seqs = 'tr_pep';
    p1_seqs = 'tr_pep';
    p2_seqs = 'tr_pep';
    
    seq_div_html = s.innerHTML;
    
    tool = new tool_pencil();
    
    c.addEventListener('mousedown', ev_canvas, false);
    c.addEventListener('mousemove', ev_canvas, false);
    c.addEventListener('mouseup',   ev_canvas, false);
}
END

    # add main jscript stuff
    $html .= $self->jscript_html();
    
    # add transcript render jscript
    $html .= $self->render_transcript_html($tr, $canvas_width, $canvas_height);

    # close javascript
    $html .= "\n</script>\n";

    $html .= "\n</body>\n</html>";
    
    open OUT, ">".$tr->stable_id.".html";
    print OUT $html;
    close OUT;
    
    # write transcript data to cache
    my $gene = $tr->{_gene};
    
    $self->{genes}->{$gene->stable_id} ||= {
        coords      => {
            chr => $tr->slice->seq_region_name,
            start => $gene->start,
            end => $gene->end,
        },
        hgnc        => $tr->{_gene_hgnc},
        transcripts => [],
    };
    push @{$self->{genes}->{$gene->stable_id}->{transcripts}}, {
        transcript => $tr->stable_id,
        mismatches => $worst_mismatches,
        length_ratio => $worst_length_ratio,
        worst_con => $worst_con->label,
    };
    
    return {};
}

sub max {
    return (sort {$a <=> $b} @_)[-1];
}

sub min {
    return (sort {$a <=> $b} @_)[0];
}

sub print_aligned {
    my $self = shift;
    my ($s1, $s2) = @_;
    
    my $width = 60;
    
    while(1) {
        my $c1 = substr($s1, 0, $width);
        my $c2 = substr($s2, 0, $width);
        
        #for my $i(0..(CORE::length($c1) - 1)) {
        #   substr($c2, $i, 1) = "." if substr($c2, $i, 1) eq substr($c1, $i, 1); 
        #}
        
        #$c1 =~ s/([A-Z]{3})/$1 /g;
        #$c2 =~ s/([A-Z]{3})/$1 /g;
        
        print "$c1\n$c2\n\n";
        last if CORE::length($s1) <= $width;
        
        $s1 = substr($s1, $width);
        $s2 = substr($s2, $width);
        
        last unless $s1 && $s2;
    }
}

sub get_ref_transcript_features {
    my $self = shift;
    my $tr   = shift;
    
    my $exons = $tr->get_all_Exons();
    return {} unless defined $exons && scalar @$exons;
    
    my $mapper = $tr->{_variation_effect_feature_cache}->{mapper};
    return {} unless defined $mapper;
    
    my $translateable_seq = $tr->{_variation_effect_feature_cache}->{translateable_seq};
    
    #$DBsingle = 1;
    
    my @features = ();
    my $tr_seq   = "";
    
    my $prev_exon;
    foreach my $exon(@$exons) {
        
        # intron
        if(defined($prev_exon)) {
            my $seq = "_" x ($exon->start > $prev_exon->start ? $exon->start - $prev_exon->end : $prev_exon->start - $exon->end);
            my $feat = {
                t => 'intron',
                s => length($tr_seq),
                seq => $seq
            };
            $tr_seq .= $seq;
            $feat->{e} = length($tr_seq);
            push @features, $feat;
        }
        $prev_exon = $exon;
        
        # map exon to CDS coords
        my @coords = $mapper->genomic2cds($exon->start, $exon->end, $exon->strand);
        
        foreach my $coord(@coords) {
            # UTR
            if($coord->isa('Bio::EnsEMBL::Mapper::Gap')) {
                my $seq = '#' x $coord->length;
                my $feat = {
                    t => 'utr',
                    s => length($tr_seq),
                    seq => $seq
                };
                $tr_seq .= $seq;
                $feat->{e} = length($tr_seq);
                push @features, $feat;
            }
            else {
                my $seq = substr($translateable_seq, $coord->start - 1, $coord->length);
                my $feat = {
                    t => 'coding',
                    s => length($tr_seq),
                    seq => $seq,
                };
                $tr_seq .= $seq;
                $feat->{e} = length($tr_seq);
                push @features, $feat;
            }
        }
    }
    
    return \@features;
}

sub get_individual_features {
    my $self         = shift;
    my $tr           = shift;
    my $ref_features = shift;
    
    # construct sequence from ref features
    my $ref_tr_seq = "";
    $ref_tr_seq .= $_->{seq} for @$ref_features;
    
    # back up features
    my @backup_features;
    
    foreach(@$ref_features) {
        push @backup_features, {
            seq => $_->{seq},
            t   => $_->{t}
        };
    }
    
    my @ind_features_array;
    
    # copy features
    foreach my $phase(0..1) {
        
        my ($new_tr_seq, $new_features);
        
        $new_tr_seq = $ref_tr_seq;
        
        # clone features
        foreach(@backup_features) {
            push @$new_features, {
                seq => $_->{seq},
                t   => $_->{t}
            };
        }
        
        push @ind_features_array, $new_features;
    }
    
    # copy for var track
    my @var_features;
    foreach(@backup_features) {
        push @var_features, {
            'seq' => $_->{t} eq 'coding' ? (' ' x length($_->{seq})) : $_->{seq},
            't' => $_->{t}
        }
    }
    
    
    foreach my $vf(sort {$b->{start} <=> $a->{start}} @{$tr->{vfs}}) {
        
        # get coords
        my ($tr_vf_start, $tr_vf_end);
        
        if($tr->strand == 1) {
            $tr_vf_start = ($vf->start - $tr->start) + 1;
            $tr_vf_end   = ($vf->end - $tr->start) + 1;
        }
        else {
            $tr_vf_start = ($tr->end - $vf->end) + 1;
            $tr_vf_end   = ($tr->end - $vf->start) + 1;
        }
        
        my $vf_length = ($vf->{end} - $vf->{start}) + 1;
        
        # get alleles
        next unless defined $vf->{genotype};
        my @alleles = @{$vf->{genotype}};
        next unless scalar @alleles;
        
        # find longest
        my $longest = (sort {length($a) <=> length($b)} @alleles)[-1];
        
        my $longest_length = length($longest);
        $longest_length = $vf_length if $vf_length > $longest_length;
        
        # pad alleles
        $_ .= "-" x ($longest_length - length($_)) for @alleles;
        
        # replace sequence
        for my $i(0..(scalar @$ref_features - 1)) {
            
            my $ref_feat = $ref_features->[$i];
            my $var_feat = $var_features[$i];
            
            if(overlap($ref_feat->{s}, $ref_feat->{e}, $tr_vf_start, $tr_vf_end)) {
             
                my ($s, $e) = ($tr_vf_start - $ref_feat->{s}, $tr_vf_end - $ref_feat->{s});
                
                # pad ref sequence if necessary
                if($longest_length > $vf_length) {
                    my $char = $ref_feat->{t} eq 'coding' ? '-' : ($ref_feat->{t} eq 'intron' ? '_' : '#');
                    substr($ref_feat->{seq}, $s, 0) = $char x ($longest_length - $vf_length);
                    substr($var_feat->{seq}, $s, 0) = ($ref_feat->{t} eq 'coding' ? 'i' : $char) x ($longest_length - $vf_length);
                }
                
                else {
                    substr($var_feat->{seq}, $s, ($e - $s) + 1) = ($vf->allele_string =~ /\-/ ? 'd' : 'v') x (($e - $s) + 1)
                }
                
                for my $phase(0..1) {
                    my $feat = $ind_features_array[$phase]->[$i];
                    my $tmp_allele = $alleles[$phase];
                    
                    my ($s, $e) = ($tr_vf_start - $ref_feat->{s}, $tr_vf_end - $ref_feat->{s});
                    
                    if($s < 0 && $tmp_allele ne '') {
                        $tmp_allele = substr($tmp_allele, 0 - $s);
                        $s = 0;
                    }
                    if($e > ($ref_feat->{e} - $ref_feat->{s})) {
                        $tmp_allele = substr($tmp_allele, 0, $e - ($ref_feat->{e} - $ref_feat->{s}));
                        $e = ($ref_feat->{e} - $ref_feat->{s});
                    }
                    
                    $tmp_allele =~ s/./\_/g if $ref_feat->{t} eq 'intron';
                    $tmp_allele =~ s/./\#/g if $ref_feat->{t} eq 'utr';
                    
                    substr($feat->{seq}, $s, ($e - $s) + 1) = $tmp_allele;
                }
            }
        }
    }
    
    push @ind_features_array, \@var_features;
    
    return \@ind_features_array;
}

sub markup_codons {
    my $self     = shift;
    my $features = shift;
    
    my $phase_adjust = 0;
    my $tmp_phase = 0;
    my $flipper = 1;
    
    for my $i(0..(scalar @$features - 1)) {
        my $feat = $features->[$i];
        $feat->{marked_seq} = $feat->{seq};
        
        if($feat->{t} eq 'coding') {
            my $j = 0;
            
            while($j < length($feat->{marked_seq})) {
                my $l = 3 - $phase_adjust;
                my $tl = $l;
                
                die("Something's gone wrong in sequence markup\n") if $l < 0;
                
                my $test = substr($feat->{marked_seq}, $j, $tl);
                my $prev_length = length($test);
                $test =~ s/\-//g;
                my $t_length = length($test);
                
                while($t_length < $l) {
                    $tl++;
                    $test = substr($feat->{marked_seq}, $j, $tl);
                    
                    last if $prev_length == length($test);
                    $prev_length = length($test);
                    $test =~ s/\-//g;
                    $t_length = length($test);
                }
                
                if($flipper) {
                    $test = substr($feat->{marked_seq}, $j, $tl);
                    $test =~ tr/ACGT-/VWXYZ/;
                    substr($feat->{marked_seq}, $j, $tl) = $test;
                }
                
                $flipper = 1 - $flipper;
                
                $j += $tl;
                
                $phase_adjust = 0;
                $tmp_phase = length($test) == 3 ? 0 : length($test);
            }
            
            $phase_adjust = $tmp_phase;
            $flipper = 1 - $flipper if $phase_adjust;
        }
    }
}

sub translate_features {
    my $self        = shift;
    my $features    = shift;
    my $codon_table = shift;
    
    my $tr_seq = "";
    $tr_seq .= $_->{seq} for @$features;
    
    $tr_seq =~ s/\#|\_|n|i|\-//g;
    
    my $pep_seq = Bio::Seq->new(-seq => $tr_seq)->translate(undef, undef, undef, $codon_table);
    
    # get sequence as 3-letter AA codes
    my $seq3 = Bio::SeqUtils->seq3($pep_seq);
    
    # shorten the sequence if a new terminator has been introduced
    while($seq3 =~ m/Ter/g) {
        my $pos = pos($seq3);
        substr($seq3, $pos) = "_" x (length($seq3) - $pos) if $pos < length($seq3);
        last;
    }
    
    return $seq3;
}

sub align_pep_to_features {
    my $self     = shift;
    my $pep      = shift;
    my $features = shift;
    
    my $tmp_pep = $pep;
    
    foreach my $f(@$features) {
        if($f->{t} eq 'coding') {
            my $tmp = $f->{seq};
            $tmp =~ s/\-//g;
            
            my $tmp2 = substr($tmp_pep, 0, length($tmp));
            
            # re-insert gaps
            for(my $i=0; $i < length($f->{seq}); $i++) {
                substr($tmp2, $i, 0) = '-' if substr($f->{seq}, $i, 1) eq '-';
            }
            
            $f->{pep} = $tmp2;
            
            # shorten tmp_pep
            $tmp_pep = substr($tmp_pep, length($tmp)) unless length($tmp_pep) < length($tmp);
        }
        
        else {
            $f->{pep} = "_" x length($f->{seq});
        }
    }
}

sub format_vfs {
    my $self = shift;
    my $tr   = shift;
    
    my @vfs;
    
    foreach my $vf(@{$tr->{vfs}}) {
        
        foreach my $tv(@{$vf->get_all_TranscriptVariations([$tr])}) {
            
            foreach my $tva(@{$tv->get_all_TranscriptVariationAlleles}) {
                push @vfs, [
                    $vf->variation_name,
                    $vf->{chr}, $vf->{start}, $vf->{end},
                    (join ",", map {$_->label} @{$tva->get_all_OverlapConsequences})
                ];
            }
        }
    }
    
    my $return = '[';
    $return .= join ",", map {"[".(join ",", map {$_ =~ /^\d+$/ ? $_ : "'$_'"} @$_)."]"} @vfs;
    $return .= ']';
    
    #print $return."\n";
    return $return;
}

sub get_karyotype_bands {
    my $self = shift;
    my $db   = shift;
    
    return if defined($self->{karyotype_bands});
    
    my $config = $self->config;
    my $ktba   = $db->get_KaryotypeBandAdaptor;
    my $sa     = $db->get_SliceAdaptor;
    
    if(!defined($ktba) || !defined($sa)) {
        warn("WARNING: VEPbrowse plugin could not get karyotype bands");
        $self->{karyotype_bands} = {};
        return;
    }
    
    foreach my $slice(grep {$_->is_reference} @{$sa->fetch_all('chromosome')}) {
        foreach my $ktb(@{$ktba->fetch_all_by_Slice($slice)}) {
            push @{$self->{karyotype_bands}->{$slice->seq_region_name}}, {
                start => $ktb->start,
                end   => $ktb->end,
                stain => $ktb->stain,
            };
        }
    }
}

sub jscript_html {
    my $self = shift;
    
    my $html =<<END;
    
function tool_pencil () {
    var tool = this;
    this.started = false;
    
    var start_x;
    
    // This is called when you start holding down the mouse button.
    // This starts the pencil drawing.
    this.mousedown = function (ev) {
        start_x = ev._x;
        tool.started = true;
    };
    
    // This function is called every time you move the mouse. Obviously, it only 
    // draws if the tool.started state is set to true (when you are holding down 
    // the mouse button).
    this.mousemove = function (ev) {
        if (tool.started) {
            var x = start_x < ev._x ? start_x : ev._x;
            var w = ev._x - start_x;
            if(w < 0) {
                w = 0 - w;
            }
            
            select_from = start_x < ev._x ? start_x : ev._x;
            select_to   = start_x > ev._x ? start_x : ev._x;
            
            // clear canvas
            ctx.clearRect(0, 0, canvas_width, canvas_height);
            
            // redraw transcript
            drawTranscript();
            
            draw_box(x, w);
        }
    };
    
    // This is called when you release the mouse button.
    this.mouseup = function (ev) {
        if (tool.started) {
            tool.mousemove(ev);
            tool.started = false;
            
            render_seq();
        }
    }
}

function ev_canvas (ev) {

    var element = c, offsetX = 0, offsetY = 0, div = document.getElementById('canvas_holder');
 
    // Compute the total offset
    if (element.offsetParent !== undefined) {
        do {
            offsetX += element.offsetLeft;
            offsetY += element.offsetTop;
        } while ((element = element.offsetParent));
    }
    
    // Add padding and border style widths to offset
    // Also add the <html> offsets in case there's a position:fixed bar
    //offsetX += div.stylePaddingLeft + div.styleBorderLeft + div.htmlLeft;
    //offsetY += div.stylePaddingTop + div.styleBorderTop + div.htmlTop;
    
    ev._x = ev.pageX - offsetX;
    ev._y = ev.pageY - offsetY;
    
    // old way
    //if (ev.layerX || ev.layerX == 0) { // Firefox
    //    ev._x = ev.layerX;
    //    ev._y = ev.layerY;
    //} else if (ev.offsetX || ev.offsetX == 0) { // Opera
    //    ev._x = ev.offsetX;
    //    ev._y = ev.offsetY;
    //}
    
    // Call the event handler of the tool.
    var func = tool[ev.type];
    if (func) {
        func(ev);
    }
}

function debug_msg(msg) {
    debug_div.innerHTML = debug_div.innerHTML + "<br/>" + msg;
    debug_div.scrollTop = debug_div.scrollHeight;
}

function toggle_value(v) {
    v = 1 - v;
    return v;
}

function draw_box(x, w) {
    
    // reset styles
    ctx.strokeStyle = "#9933CC";
    ctx.fillStyle = "rgba(153, 51, 204, 0.2)";
    
    // draw box
    ctx.fillRect(x, 0, w, canvas_height);
    ctx.strokeRect(x, 0, w, canvas_height);
}

function select_all() {
    
    select_from = 0;
    select_to   = canvas_width;
    
    // clear canvas
    ctx.clearRect(0, 0, canvas_width, canvas_height);
    
    // redraw transcript
    drawTranscript();
    
    draw_box(0, canvas_width);
    
    render_seq();
}

function clear_selection() {
    
    select_from = 0;
    select_to   = 0;
    
    // clear canvas
    ctx.clearRect(0, 0, canvas_width, canvas_height);
    
    // redraw transcript
    drawTranscript();
    
    render_seq();
}

function zoom(factor) {
    var mid_point = (select_to + select_from) / 2;
    var dist = mid_point - select_from;
    
    dist = dist * factor;
    
    select_from = mid_point - dist;
    select_to = mid_point + dist;
    
    // clear canvas
    ctx.clearRect(0, 0, canvas_width, canvas_height);
    
    // redraw transcript
    drawTranscript();
    
    draw_box(select_from, (2 * dist));
    
    render_seq();
}

function move(factor) {
    dist = (select_to - select_from + 1);
    var amount = dist * factor;
    
    select_from = select_from + amount;
    select_to = select_to + amount;
    
    // clear canvas
    ctx.clearRect(0, 0, canvas_width, canvas_height);
    
    // redraw transcript
    drawTranscript();
    
    draw_box(select_from, dist);
    
    render_seq();
}
END

    return $self->render_seq_jscript.$html;
}

sub render_seq_jscript {
    my $self = shift;
    
    my $html =<<END;
    
function render_seq () {
    var i=0;
    var rendered = '';
    
    // convert coords to bp
    var bp_from, bp_to;
    
    bp_from = Math.round((select_from - x_off) / x_scale);
    bp_to = Math.round((select_to - x_off) / x_scale);
    
    if(bp_from < 1) {
        bp_from = 1;
    }
    if(bp_to < 1) {
        bp_to = 1;
    }
    
    if(bp_from == bp_to) {
        s.innerHTML = seq_div_html;
        return;
    }
    
    debug_msg("Rendering FROM " + bp_from + " TO " + bp_to);
    
    // fix for gaps in ref align
    for(i=0; i<gaps.length; i++) {
        if(bp_from > gaps[i]) {
            bp_from++;
            debug_msg("Added 1 to bp_from");
        }
        if(bp_to > gaps[i]) {
            bp_to++;
            debug_msg("Added 1 to bp_to");
        }
    }
    
    rendered = rendered + "<p>FROM " + bp_from + " TO " + bp_to + "</p>";
    
    // clone seqs array
    var tmp_seqs = seqs.slice(0);
    
    for(i=0; i<seqs.length; i++) {
        tmp_seqs[i] = tmp_seqs[i].substr(bp_from, (bp_to - bp_from) + 1);
    }
    
    // remove non-coding parts if requested
    if(coding_only) {
        for(i=0; i<tmp_seqs.length; i++) {
            tmp_seqs[i] = tmp_seqs[i].replace(/[\_\#]/g, "");
        }
    }
    
    for(i=0; i<=tmp_seqs[0].length;i+=60) {
        var coord = bp_from + i;
        
        var pre;
        pre = coord + " ";
        
        var j;
        for(j=String(coord).length; j<=String(seqs[0].length).length; j++) {
            pre = pre + " ";
        }
        
        var line_seq = tmp_seqs[0].substr(i, 60);
        var after_coord = coord + (line_seq.length - 1);
        
        var post = " ";
        
        for(j=line_seq.length; j<=60; j++) {
            post = post + " ";
        }
        
        // render sequences
        for(var j=0; j<tmp_seqs.length; j++) {
            var seq = tmp_seqs[j].substr(i, 60);
            
            var is_pep = (j % 2 == 0 ? 0 : 1);
            var label, coord1, coord2, seq_name;
            coord1 = "";
            coord2 = "";
            
            seq_name = (j < 2 ? "REF" : "PH" + Math.floor(j / 2));
            
            if(j == (tmp_seqs.length - 1)) {
                seq_name = "VAR";
            }
            
            if(!is_pep && j != 0) {
                rendered = rendered + "<hr class='seq_spacer'/>";
            }
            
            // are we rendering this sequence?
            if(seq_name == "REF") {
                if(is_pep) {
                    if(!ref_seqs.match(/pep/)) {
                        continue;
                    }
                }
                else {
                    if(!ref_seqs.match(/tr/)) {
                        continue;
                    }
                }
            }
            else if(Math.floor(j / 2) == 1) {
                if(is_pep) {
                    if(!p1_seqs.match(/pep/)) {
                        continue;
                    }
                }
                else {
                    if(!p1_seqs.match(/tr/)) {
                        continue;
                    }
                }
            }
            else if(Math.floor(j / 2) == 2) {
                if(is_pep) {
                    if(!p2_seqs.match(/pep/)) {
                        continue;
                    }
                }
                else {
                    if(!p2_seqs.match(/tr/)) {
                        continue;
                    }
                }
            }
            
            if(j == (tmp_seqs.length - 1)) {
                label = seq_name + "     | ";
                for(var k=0; k<pre.length; k++) {
                    coord1 = coord1 + " ";
                }                
                coord2 = "";
                seq = format_var_seq(seq, bp_from);
            }
            else if(is_pep) {
                label = seq_name + " PEP | ";
                for(var k=0; k<pre.length; k++) {
                    coord1 = coord1 + " ";
                }                
                coord2 = "";
                seq = format_pep_seq(seq);
            }
            else {
                label = seq_name + " TR  | ";
                coord1 = pre;
                coord2 = after_coord;
                seq = format_tr_seq(seq);
            }
            
            rendered = rendered + label + coord1 + seq + post + coord2 + "<br/>";
        }
        
        rendered = rendered + "<br/><br/>";
    }
    
    s.innerHTML = '<pre>'+rendered+'</pre>';
}

function format_tr_seq(seq) {
    seq = seq.replace(/\_/g, '<span class="intron" title="Intron">i</span>');
    seq = seq.replace(/A/g, '<span class="base_a codon_1">A</span>');
    seq = seq.replace(/C/g, '<span class="base_c codon_1">C</span>');
    seq = seq.replace(/G/g, '<span class="base_g codon_1">G</span>');
    seq = seq.replace(/T/g, '<span class="base_t codon_1">T</span>');
    seq = seq.replace(/-/g, '<span class="codon_1">-</span>');
    seq = seq.replace(/V/g, '<span class="base_a codon_2">A</span>');
    seq = seq.replace(/W/g, '<span class="base_c codon_2">C</span>');
    seq = seq.replace(/X/g, '<span class="base_g codon_2">G</span>');
    seq = seq.replace(/Y/g, '<span class="base_t codon_2">T</span>');
    seq = seq.replace(/Z/g, '<span class="codon_2">-</span>');
    seq = seq.replace(/\#/g, '<span class="utr" title="UTR">u</span>');

    return seq;
}

function format_pep_seq(seq) {
    seq = seq.replace(/\#/g, '<span class="utr">u</span>');
    seq = seq.replace(/\_/g, ' ');
    seq = seq.replace(/(M\-*e\-*t)/g, '<span class="start_codon">\$1</span>');
    seq = seq.replace(/(T\-*e\-*r)/g, '<span class="stop_codon">\$1</span>');
    
    seq = seq.replace(/(G\-*l\-*y)/g, '<span class="pep_1">\$1</span>');
    seq = seq.replace(/(P\-*r\-*o)/g, '<span class="pep_1">\$1</span>');
    seq = seq.replace(/(S\-*e\-*r)/g, '<span class="pep_1">\$1</span>');
    seq = seq.replace(/(T\-*h\-*r)/g, '<span class="pep_1">\$1</span>');
    
    seq = seq.replace(/(H\-*i\-*s)/g, '<span class="pep_2">\$1</span>');
    seq = seq.replace(/(L\-*y\-*s)/g, '<span class="pep_2">\$1</span>');
    seq = seq.replace(/(A\-*r\-*g)/g, '<span class="pep_2">\$1</span>');
    
    seq = seq.replace(/(P\-*h\-*e)/g, '<span class="pep_3">\$1</span>');
    seq = seq.replace(/(T\-*r\-*p)/g, '<span class="pep_3">\$1</span>');
    seq = seq.replace(/(T\-*y\-*r)/g, '<span class="pep_3">\$1</span>');
    
    seq = seq.replace(/(I\-*l\-*e)/g, '<span class="pep_4">\$1</span>');
    seq = seq.replace(/(L\-*e\-*u)/g, '<span class="pep_4">\$1</span>');
    seq = seq.replace(/(M\-*e\-*t)/g, '<span class="pep_4">\$1</span>');
    seq = seq.replace(/(V\-*a\-*l)/g, '<span class="pep_4">\$1</span>');

    return seq;
}

function format_var_seq(seq, bp_from) {
    for(var i=0; i<vfs.length; i++) {
        var name = vfs[i][0], chr = vfs[i][1], start = vfs[i][2], end = vfs[i][3], con = vfs[i][4];
        
        //continue if start < bp_from;
    }
    
    seq = seq.replace(/\#/g, ' ');
    seq = seq.replace(/\_/g, ' ');

    return seq;
}
END
    return $html;
}

sub render_transcript_html {
    my $self          = shift;
    my $tr            = shift;
    my $canvas_width  = shift;
    my $canvas_height = shift;
    
    my $exons = $tr->get_all_Exons();
    
    my $tr_start = $tr->start;
    my $tr_end   = $tr->end;
    my $width    = $tr_end - $tr_start + 1;
    
    my $x_scale  = ($canvas_width - 20)  / $width;
    my $y_scale  = ($canvas_height - 30) / 100;
    my $x_off    = 10;
    my $y_off    = 10;
    
    my $zero_string = '0' x (CORE::length(int(100 / $x_scale)) - 1);
    my $bases_per_bar = '1'.$zero_string;
    
    my $start = int($tr_start / $bases_per_bar) * $bases_per_bar;
    my $end = (int($tr_start / $bases_per_bar) + 1) * $bases_per_bar;
    my $colour = 0;
    
    my $html = "";
    
    $html .= "\n\n// scale\n";
    $html .= "ctx.strokeStyle=\"#c3c3c3\";\n";
    $html .= "ctx.font=\"8px Arial\";\n";
    $html .= "ctx.beginPath();\n";
    
    my ($x, $w);
    
    while($start < $tr_end) {
        my $method = $colour ? 'fillRect' : 'strokeRect';
        
        $x = $x_off + (($start - $tr_start) * $x_scale);
        $w = (($end - $start) * $x_scale); 
        $x = $canvas_width - ($x + $w) if $tr->strand == -1;
        
        $html .= "ctx.fillStyle=\"#c3c3c3\";\n";
        $html .= "ctx.$method(".join(",", (
            $x,
            $canvas_height - 15,
            $w,
            5,
        )).");\n";
        
        # tick and label
        if($start =~ /(5|0)$zero_string$/) {
            
            # strand indicator
            my $fill = $colour ? '#FFFFFF' : '#C3C3C3';
            
            $x = $x + $w if $tr->strand == -1;
            
            $html .= "ctx.beginPath();\n";
            $html .= "ctx.moveTo(".join(",", (
                $x,
                $canvas_height - 15,
            )).");\n";
            
            if($tr->strand == 1) {
                $html .= "ctx.lineTo(".join(",", (
                    $x + 10,
                    $canvas_height - 12.5,
                )).");\n";
            }
            else {
                $html .= "ctx.lineTo(".join(",", (
                    $x - 10,
                    $canvas_height - 12.5,
                )).");\n";
            }
            
            $html .= "ctx.lineTo(".join(",", (
                $x,
                $canvas_height - 10,
            )).");\n";
            $html .= "ctx.lineTo(".join(",", (
                $x,
                $canvas_height - 15,
            )).");\n";
            $html .= "ctx.stroke();\n";
            $html .= "ctx.fillStyle=\"$fill\";\n";
            $html .= "ctx.fill();\n";
            $html .= "ctx.closePath();\n";
            
            # coords and guide lines
            my $string = $start;
            1 while $string =~ s/^(-?\d+)(\d\d\d)/$1,$2/;
            
            $x = $x_off + (($start - $tr_start) * $x_scale);
            $x = $canvas_width - $x if $tr->strand == -1;
            
            $html .= "ctx.fillStyle=\"#000000\";\n";
            $html .= "ctx.fillText(\"$string\",".join(",", (
                $x + 2,
                $canvas_height - 20
            )).");\n";
            
            $html .= "ctx.strokeStyle=\"#E8E8E8\";\n";
            $html .= "ctx.beginPath();\n";
            $html .= "ctx.moveTo(".join(",", (
                $x,
                0,
            )).");\n";
            $html .= "ctx.lineTo(".join(",", (
                $x,
                $canvas_height,
            )).");\n";
            $html .= "ctx.stroke();\n";
            $html .= "ctx.closePath();\n";
        }
        
        $colour = 1 - $colour;
        $start += $bases_per_bar;
        $end += $bases_per_bar;
    }
    
    $html .= "ctx.closePath();\n";
    
    # render introns
    $html .= "\n\n// introns\n";
    $html .= "ctx.beginPath();\n";
    $html .= "ctx.strokeStyle=\"#0099FF\";\n";
    
    foreach my $intron(@{$tr->get_all_Introns}) {
        
        $x = $x_off + (($intron->start - $tr_start) * $x_scale);
        $x = $canvas_width - $x if $tr->strand == -1;
        $html .= "ctx.moveTo(".join(",", (
            $x,
            $y_off + (20 * $y_scale)
        )).");\n";
        
        $x = $x_off + (((($intron->start + $intron->end) / 2) - $tr_start) * $x_scale);
        $x = $canvas_width - $x if $tr->strand == -1;
        $html .= "ctx.lineTo(".join(",", (
            $x,
            $y_off + (10 * $y_scale)
        )).");\n";
        
        $x = $x_off + (($intron->end - $tr_start) * $x_scale);
        $x = $canvas_width - $x if $tr->strand == -1;
        $html .= "ctx.lineTo(".join(",", (
            $x,
            $y_off + (20 * $y_scale)
        )).");\n";
        
        $html .= "ctx.stroke();\n";
    }
    
    $html .= "\n\n// exons\n";
    $html .= "ctx.closePath();\n";
    $html .= "ctx.fillStyle=\"#0000FF\";\n";
    
    # render exons
    foreach my $exon(@$exons) {
        
        # non-coding part
        $x = $x_off + (($exon->start - $tr_start) * $x_scale);
        $w = (($exon->end - $exon->start) * $x_scale); 
        $x = $canvas_width - ($x + $w) if $tr->strand == -1;
        
        $html .= "ctx.strokeRect(".join(",", (
            $x,
            $y_off + (10 * $y_scale),
            $w,
            20 * $y_scale,
        )).");\n";
        
        # coding part
        my ($cs, $ce) = ($exon->coding_region_start($tr), $exon->coding_region_end($tr));
        
        next unless (defined($cs) and defined($ce));
        
        $x = $x_off + (($cs - $tr_start) * $x_scale);
        $w = (($ce - $cs) * $x_scale); 
        $x = $canvas_width - ($x + $w) if $tr->strand == -1;
        
        $html .= "ctx.fillRect(".join(",", (
            $x - 1,
            $y_off + (0 * $y_scale),
            $w + 2,
            40 * $y_scale,
        )).");\n";
    }
    
    # render variants
    $html .= "\n\n// variants\n";
    
    foreach my $vf(@{$tr->{vfs}}) {
        my $ref_allele = (split /\//, $vf->allele_string)[0];
        
        my $y = $y_off + (60 * $y_scale);
        
        foreach my $allele(@{$vf->{genotype}}) {
            my $colour = $allele eq $ref_allele ? '00FF00' : 'FF0000';
            
            $x = $x_off + (($vf->start - $tr_start) * $x_scale);
            $w = (($vf->end - $vf->start) * $x_scale);
            $x = $canvas_width - ($x + $w) if $tr->strand == -1;
            
            if($w < 0) {
                $x += $w;
                $w = 0 - $w;
            }
            $w = 1 if $w < 1;
            
            $html .= "ctx.fillStyle=\"#$colour\";\n";
            $html .= "ctx.fillRect(".join(",", (
                $x,
                $y,
                $w,
                10 * $y_scale,
            )).");\n";
            
            $y += (10 * $y_scale);
        }
    }
    
    $html =~ s/(^|\n)/$1\t/g;
    $html = "function drawTranscript(from, to) {\n$html\n}\n";
    
    return $html;
}

sub render_karyotype_html {
    my $self = shift;
    my $canvas_width = shift;
    my $canvas_height = shift;
    
    my $bands = $self->{karyotype_bands};
    my $x_off = 10;
    my $y_off = 10;
    
    $DB::single = 1;
    
    my $longest_chr = (sort {$a <=> $b} map {$_->{end}} map {@{$_}} values %{$bands})[-1];
    my $y_scale = ($canvas_height - (4 * $y_off)) / $longest_chr;
    my $chr_width = ($canvas_width - (2 * $x_off)) / ((3 * scalar keys %$bands) - 1);
    
    my $html = "\nfunction drawKaryotype() {\ngenes = new Array();";
    
    my $x = $x_off;
    
    my $stain_colours = {
        'gneg' => 'white',
        'gpos25' => 'lightgrey',
        'gpos50' => 'grey',
        'gpos75' => 'darkgrey',
        'gpos100' => 'black',
        'acen' => 'yellow',
        'gvar' => 'lightgrey',
        'stalk' => 'lightblue',
    };
    
    my %x_offs = ();
    
    foreach my $chr(sort {($a !~ /^\d+$/ || $b !~ /^\d+/) ? $a cmp $b : $a <=> $b} keys %$bands) {
        my $length = $bands->{$chr}->[-1]->{end};
        my $y_adj  = $canvas_height - ((4 * $y_off) + ($length * $y_scale));
        
        # chr outline
        $html .= "ctx.strokeStyle=\"black\";\n";
        $html .= "ctx.strokeRect(".join(",", (
            $x,
            $y_off + $y_adj,
            $chr_width,
            $length * $y_scale,
        )).");\n";
        
        #$html .= "ctx.restore();\n";
        
        # chr label
        $html .= "ctx.fillStyle=\"#000000\";\n";
        $html .= "ctx.fillText(\"$chr\",".join(",", (
            $x + 2,
            $canvas_height - 10
        )).");\n";
        
        
        # karyotype bands
        my $colour = 'white';
        
        foreach my $band(@{$bands->{$chr}}) {
            my $y = $y_off + ($band->{start} * $y_scale) + $y_adj;
            my $h = ($band->{end} - $band->{start}) * $y_scale;
            
            $colour = $stain_colours->{$band->{stain}} || 'white';
            
            $html .= "ctx.fillStyle=\"$colour\";\n";
            $html .= "ctx.fillRect(".join(",", (
                $x + 1,
                $y,
                $chr_width - 2,
                $h,
            )).");\n";
        }
        
        $x_offs{$chr} = $x + $chr_width;
        
        $x += (3 * $chr_width);
    }
    
    foreach my $gene(keys %{$self->{genes}}) {
        
        my $coords = $self->{genes}->{$gene}->{coords};
        
        my $length = $bands->{$coords->{chr}}->[-1]->{end};
        my $y_adj  = $canvas_height - ((4 * $y_off) + ($length * $y_scale));
        
        my $x = $x_offs{$coords->{chr}};
        my $y = $y_off + (((($coords->{end} + $coords->{start}) / 2) * $y_scale)) + $y_adj;
        
        my $trs = $self->{genes}->{$gene}->{transcripts};
        my ($lr, $mm) = (100, 0);
        
        foreach my $tr(@$trs) {
            $lr = $tr->{length_ratio} if $tr->{length_ratio} < $lr;
            $mm = $tr->{mismatches} if $tr->{mismatches} > $mm;
        }
        
        my $colour = 'grey';
        
        # no length changes
        if($lr == 100) {
            if($mm <= 1) {
                $colour = 'green';
            }
            elsif($mm < 5) {
                $colour = 'orange';
            }
            else {
                $colour = 'red';
            }
        }
        else {
            if($lr < 50) {
                $colour = 'red';
            }
            elsif($lr < 75) {
                $colour = 'orange';
            }
            else {
                $colour = 'green';
            }
        }
        
        $html .= "ctx.fillStyle=\"$colour\";\n";
        $html .= "ctx.fillRect(".join(",", (
            $x + 10,
            $y - 3,
            10,
            6,
        )).");\n";
        
        $html .= "ctx.beginPath();\n";
        $html .= "ctx.strokeStyle=\"$colour\";\n";
        $html .= "ctx.moveTo(".join(",", (
            $x + 10,
            $y - 2
        )).");\n";
        $html .= "ctx.lineTo(".join(",", (
            $x + 2,
            $y,
        )).");\n";
        $html .= "ctx.lineTo(".join(",", (
            $x + 10,
            $y + 2
        )).");\n";
        
        $html .= "ctx.stroke();\n";
        $html .= "ctx.fill();\n";
        $html .= "ctx.closePath();\n";
        
        my $gene_data = '['.(join ",", map {"['".$_->{transcript}."',".$_->{mismatches}.",".$_->{length_ratio}.",'".$_->{worst_con}."']"} @{$trs}).']';
        
        $html .= "genes.push({".
            "name: '".$self->{genes}->{$gene}->{hgnc}."', ".
            "selected: 0, ".
            "trs: $gene_data, ".
            "g: '$gene', ".
            "c: '".($coords->{chr}.":".$coords->{start}."-".$coords->{end})."', x1: ".($x+10).", y1: ".($y-2).", x2: ".($x+20).", y2: ".($y+3)."});";
    }
    
    $html .= "}\n";
    return $html;
}

sub karyotype_jscript {
    my $self = shift;
    
    my $html =<<END;


function tool_pencil () {
    var tool = this;
    
    var x, y;
    
    this.mousemove = function (ev) {
        x = ev._x;
        y = ev._y;
        
        var overlap = 0;
        
        for(var i=0; i<genes.length; i++) {
            if(x >= genes[i].x1 && x <= genes[i].x2 && y >= genes[i].y1 && y <= genes[i].y2) {
                overlap = 1;
                break;
            }
        }
        
        if(overlap == 1) {
            c.style.cursor = 'pointer';
        }
        else {
            c.style.cursor = 'auto';
        }
    }
    
    this.mousedown = function (ev) {
        x = ev._x;
        y = ev._y;
        
        var html = "";
        
        // clear canvas
        ctx.clearRect(0, 0, canvas_width, canvas_height);
        
        drawKaryotype(canvas_width, canvas_height);
        
        for(var i=0; i<genes.length; i++) {
            if(x >= genes[i].x1 && x <= genes[i].x2 && y >= genes[i].y1 && y <= genes[i].y2) {
                genes[i].selected = 1 - genes[i].selected;
                
                if(genes[i].selected == 1) {
                    var gene = genes[i], trs = genes[i].trs;
                    
                    html = html + '<h3>' + gene.name + ' (' + gene.g + ') | ' + gene.c + '</h3>';
                    html = html + '<ul>';
                    
                    for(var j=0; j<trs.length; j++) {
                        html = html + '<li><a href="' + trs[j][0] + '.html">' + trs[j][0] + '</a>';
                        html = html + ' | Mismatches: ' + trs[j][1] + ' | Length ratio: ' + trs[j][2];
                        html = html + ' | Worst consequence: ' + trs[j][3];
                        html = html + '</li>';
                    }
                    
                    html = html + '</ul>';
                    
                    ctx.strokeStyle = 'blue';
                    ctx.strokeRect(gene.x1-1, gene.y1-1, (gene.x2 - gene.x1)+1, (gene.y2 - gene.y1)+1);
                    ctx.strokeRect(gene.x1-2, gene.y1-2, (gene.x2 - gene.x1)+3, (gene.y2 - gene.y1)+3);
                }
            }
        }
        
        info.innerHTML = html;
    };
}

function ev_canvas (ev) {

    var element = c, offsetX = 0, offsetY = 0, div = document.getElementById('canvas_holder');
 
    // Compute the total offset
    if (element.offsetParent !== undefined) {
        do {
            offsetX += element.offsetLeft;
            offsetY += element.offsetTop;
        } while ((element = element.offsetParent));
    }
    
    // Add padding and border style widths to offset
    // Also add the <html> offsets in case there's a position:fixed bar
    //offsetX += div.stylePaddingLeft + div.styleBorderLeft + div.htmlLeft;
    //offsetY += div.stylePaddingTop + div.styleBorderTop + div.htmlTop;
    
    ev._x = ev.pageX - offsetX;
    ev._y = ev.pageY - offsetY;
    
    // Call the event handler of the tool.
    var func = tool[ev.type];
    if (func) {
        func(ev);
    }
}
END
    return $html;
}

sub get_css {
    my $self = shift;
    
    my $html =<<END;
    
<style type="text/css">
    body {font-family: arial, sans-serif; }
    div {padding: 5px;}
    
    a {color: blue; text-decoration: none;}
    a.visited {color: blue; text-decoration: none;}
    

    .base_a {color:blue;}
    .base_c {color:red;}
    .base_g {color:green;}
    .base_t {color:orange;}
    
    .codon_1 {background: #fff9af;}
    .codon_2 {border-color: #ffffff;}
    
    .start_codon {background: green; color: white;}
    .stop_codon {background: red; color: white;}
    
    .pep_1 { background: orange; }
    .pep_2 { background: red; color: white; }
    .pep_3 { background: blue; color: white; }
    .pep_4 { background: green; color: white; }
    
    .seq_spacer { border: 1px dashed; color: lightgrey; margin: 2px; }
    .control_spacer { border: 1px dashed; color: black; margin: 2px; }
    
    .intron {background: lightgrey; color: grey;}
    .utr    {background: grey; color: lightgrey;}
    
    .control {
        font-family: arial, sans-serif;
        font-size: small;
        border:1px solid #c3c3c3;
        background: lightgrey;
    }
    .header {
        border:1px solid #c3c3c3;
    }
</style>

END

    return $html;
}

1;

