# -*- coding: utf-8 -*-
# Generated by Django 1.9.7 on 2016-09-30 15:18
from __future__ import unicode_literals

from django.db import migrations
import randomslugfield.fields


class Migration(migrations.Migration):

    dependencies = [
        ('base', '0005_alignment_populate_slug_values'),
    ]

    operations = [
        migrations.AlterField(
            model_name='alignment',
            name='slug',
            field=randomslugfield.fields.RandomSlugField(blank=True, editable=False,
                                                         length=16, max_length=16, unique=True),
        ),
    ]