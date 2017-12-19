# LaTeX2HTML 2012 (1.2)
# Associate labels original text with physical files.


$key = q/SubSec:para/;
$external_labels{$key} = "$URL/" . q|node7.html|; 
$noresave{$key} = "$nosave";

$key = q/SubSec:Examples/;
$external_labels{$key} = "$URL/" . q|node5.html|; 
$noresave{$key} = "$nosave";

$key = q/Sec:para/;
$external_labels{$key} = "$URL/" . q|node6.html|; 
$noresave{$key} = "$nosave";

1;


# LaTeX2HTML 2012 (1.2)
# labels from external_latex_labels array.


$key = q/SubSec:para/;
$external_latex_labels{$key} = q|4.1|; 
$noresave{$key} = "$nosave";

$key = q/SubSec:Examples/;
$external_latex_labels{$key} = q|3.1|; 
$noresave{$key} = "$nosave";

$key = q/Sec:para/;
$external_latex_labels{$key} = q|4|; 
$noresave{$key} = "$nosave";

1;

