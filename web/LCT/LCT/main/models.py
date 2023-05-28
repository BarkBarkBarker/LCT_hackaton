from django.db import models


class Material(models.Model):
    turning = "TR"
    triple = "TP"
    quadruple = "QP"
    reduction = "RQ"
    tube = "TB"
    TYPE = [
        (turning, "Отвод"),
        (triple, "Тройник"),
        (quadruple, "Крестовина"),
        (reduction, "Редукция"),
        (tube, "Труба"),
    ]
    type = models.CharField(choices=TYPE, max_length=60)
    id = models.IntegerField(primary_key=True)
    d_1 = models.IntegerField()
    d_2 = models.IntegerField()
    d_3 = models.IntegerField()
    angle = models.IntegerField()

