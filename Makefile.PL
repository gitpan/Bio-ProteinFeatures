# $Id: Makefile.PL,v 1.1.1.1 2003/09/06 06:36:16 cvspub Exp $
use ExtUtils::MakeMaker;


WriteMakefile(
    'NAME'		=> 'Bio::ProteinFeatures',
    'VERSION_FROM'	=> 'ProteinFeatures.pm', 
    'PREREQ_PM'		=> {}, 
    ($] >= 5.005 ?    
      (ABSTRACT_FROM => 'ProteinFeatures.pm', 
       AUTHOR     => 'xern <xern@cpan.org>') : ()),
);
