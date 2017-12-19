# LaTeX2HTML 2012 (1.2)
# Associate internals original text with physical files.


$key = q/SubSec:para/;
$ref_files{$key} = "$dir".q|node7.html|; 
$noresave{$key} = "$nosave";

$key = q/SubSec:Examples/;
$ref_files{$key} = "$dir".q|node5.html|; 
$noresave{$key} = "$nosave";

$key = q/Sec:para/;
$ref_files{$key} = "$dir".q|node6.html|; 
$noresave{$key} = "$nosave";

1;

